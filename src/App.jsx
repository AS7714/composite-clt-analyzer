import { useState, useMemo, useCallback } from "react";

// ============================================================
// MICROMECHANICS: Fiber + Matrix → Ply Properties
// ============================================================

function micromechanics(fiber, matrix, Vf) {
  const Vm = 1 - Vf;
  const Ef1 = fiber.Ef1, Ef2 = fiber.Ef2, Gf12 = fiber.Gf12, nuf12 = fiber.nuf12;
  const Em = matrix.Em, Gm = matrix.Gm, num = matrix.num;

  // Rule of Mixtures (Voigt) — longitudinal
  const E1 = Vf * Ef1 + Vm * Em;

  // Halpin-Tsai — transverse modulus
  const xi_E2 = 2.0; // reinforcement factor for circular fibers
  const eta_E2 = (Ef2 / Em - 1) / (Ef2 / Em + xi_E2);
  const E2 = Em * (1 + xi_E2 * eta_E2 * Vf) / (1 - eta_E2 * Vf);

  // Halpin-Tsai — shear modulus
  const xi_G = 1.0;
  const eta_G = (Gf12 / Gm - 1) / (Gf12 / Gm + xi_G);
  const G12 = Gm * (1 + xi_G * eta_G * Vf) / (1 - eta_G * Vf);

  // Rule of Mixtures — Poisson's ratio
  const nu12 = Vf * nuf12 + Vm * num;

  // Inverse ROM estimates (for comparison display)
  const E2_rom = 1 / (Vf / Ef2 + Vm / Em);
  const G12_rom = 1 / (Vf / Gf12 + Vm / Gm);

  return { E1, E2, G12, nu12, nu21: nu12 * E2 / E1, E2_rom, G12_rom };
}

// ============================================================
// CLT ENGINE
// ============================================================

function computeQ(E1, E2, G12, nu12) {
  const nu21 = nu12 * E2 / E1;
  const d = 1 - nu12 * nu21;
  return { Q11: E1/d, Q22: E2/d, Q12: nu12*E2/d, Q66: G12, nu21 };
}

function computeQbar(Q, thetaDeg) {
  const t = thetaDeg * Math.PI / 180;
  const { Q11, Q22, Q12, Q66 } = Q;
  const U1 = (3*Q11 + 3*Q22 + 2*Q12 + 4*Q66)/8;
  const U2 = (Q11 - Q22)/2;
  const U3 = (Q11 + Q22 - 2*Q12 - 4*Q66)/8;
  const U4 = (Q11 + Q22 + 6*Q12 - 4*Q66)/8;
  const U5 = (Q11 + Q22 - 2*Q12 + 4*Q66)/8;
  const c2 = Math.cos(2*t), s2 = Math.sin(2*t);
  const c4 = Math.cos(4*t), s4 = Math.sin(4*t);
  return {
    Q11b: U1+U2*c2+U3*c4, Q22b: U1-U2*c2+U3*c4,
    Q12b: U4-U3*c4, Q66b: U5-U3*c4,
    Q16b: U2*s2/2+U3*s4, Q26b: U2*s2/2-U3*s4
  };
}

function computeABD(angles, plyT, Q) {
  const n = angles.length, h = n * plyT;
  const z = Array.from({length: n+1}, (_, i) => -h/2 + i*plyT);
  const A=[[0,0,0],[0,0,0],[0,0,0]];
  const B=[[0,0,0],[0,0,0],[0,0,0]];
  const D=[[0,0,0],[0,0,0],[0,0,0]];
  const qbars = [];
  for (let k = 0; k < n; k++) {
    const qb = computeQbar(Q, angles[k]);
    qbars.push(qb);
    const qM = [[qb.Q11b,qb.Q12b,qb.Q16b],[qb.Q12b,qb.Q22b,qb.Q26b],[qb.Q16b,qb.Q26b,qb.Q66b]];
    const zk=z[k+1], zk1=z[k];
    for (let i=0;i<3;i++) for (let j=0;j<3;j++) {
      A[i][j]+=qM[i][j]*(zk-zk1);
      B[i][j]+=qM[i][j]*(zk*zk-zk1*zk1)/2;
      D[i][j]+=qM[i][j]*(zk*zk*zk-zk1*zk1*zk1)/3;
    }
  }
  return { A,B,D,z,qbars,h,n };
}

function inv3(m) {
  const [[a,b,c],[d,e,f],[g,h,i]]=m;
  const det=a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g);
  if(Math.abs(det)<1e-30) return null;
  return [
    [(e*i-f*h)/det,(c*h-b*i)/det,(b*f-c*e)/det],
    [(f*g-d*i)/det,(a*i-c*g)/det,(c*d-a*f)/det],
    [(d*h-e*g)/det,(b*g-a*h)/det,(a*e-b*d)/det]
  ];
}

function effInplane(A, h) {
  const a=inv3(A); if(!a) return null;
  return { Ex:1/(h*a[0][0]), Ey:1/(h*a[1][1]), Gxy:1/(h*a[2][2]), nuxy:-a[0][1]/a[0][0] };
}
function effFlexural(D, h) {
  const d=inv3(D); if(!d) return null;
  const h3=h*h*h;
  return { Efx:12/(h3*d[0][0]), Efy:12/(h3*d[1][1]), Gfxy:12/(h3*d[2][2]), nufxy:-d[0][1]/d[0][0] };
}

function solveResponse(A,B,D,N,M) {
  // Build 6x6 and solve via simple Gaussian elimination for small system
  const S = Array.from({length:6}, (_,i) => Array.from({length:7}, (_,j) => {
    if (j===6) return i<3?N[i]:M[i-3];
    const r=i<3?i:i-3, c=j<3?j:j-3;
    if(i<3&&j<3) return A[r][c];
    if(i<3&&j>=3) return B[r][c];
    if(i>=3&&j<3) return B[r][c];
    return D[r][c];
  }));
  // Gaussian elimination
  for(let col=0;col<6;col++){
    let maxR=col;
    for(let r=col+1;r<6;r++) if(Math.abs(S[r][col])>Math.abs(S[maxR][col])) maxR=r;
    [S[col],S[maxR]]=[S[maxR],S[col]];
    if(Math.abs(S[col][col])<1e-30) continue;
    for(let r=col+1;r<6;r++){
      const f=S[r][col]/S[col][col];
      for(let c=col;c<=6;c++) S[r][c]-=f*S[col][c];
    }
  }
  const x=new Array(6).fill(0);
  for(let i=5;i>=0;i--){
    let s=S[i][6];
    for(let j=i+1;j<6;j++) s-=S[i][j]*x[j];
    x[i]=s/S[i][i];
  }
  return { eps0:x.slice(0,3), kappa:x.slice(3,6) };
}

function plyStresses(qbars, z, eps0, kappa, angles) {
  const results = [];
  for (let k=0; k<qbars.length; k++) {
    const qb=qbars[k];
    const qM=[[qb.Q11b,qb.Q12b,qb.Q16b],[qb.Q12b,qb.Q22b,qb.Q26b],[qb.Q16b,qb.Q26b,qb.Q66b]];
    const zB=z[k],zT=z[k+1],zM=(zB+zT)/2;
    const pts=[zB,zM,zT].map(zv => {
      const e=[eps0[0]+zv*kappa[0], eps0[1]+zv*kappa[1], eps0[2]+zv*kappa[2]];
      const s=[0,1,2].map(i => (qM[i][0]*e[0]+qM[i][1]*e[1]+qM[i][2]*e[2])*1000); // MPa
      return {z:zv, strain:e, stress:s};
    });
    results.push({ply:k, angle:angles[k], zB, zT, pts});
  }
  return results;
}

// ============================================================
// COLORS & STYLING
// ============================================================
const C = {
  bg:'#0a0e17', sf:'#111827', sf2:'#1a2235', bd:'#2a3650', bd2:'#3a4a65',
  tx:'#e2e8f0', dim:'#8896b0', ac:'#38bdf8', ac2:'#818cf8', pk:'#f472b6',
  gn:'#34d399', or:'#fb923c', rd:'#f87171', yl:'#fbbf24',
};
const PLY_C = {0:'#38bdf8',90:'#f472b6',45:'#34d399','-45':'#fb923c'};
function plyColor(a) {
  const n=((a%360)+360)%360;
  return PLY_C[n]||PLY_C[360-n]||PLY_C[n-360]||'#8896b0';
}
const mono = "'JetBrains Mono', monospace";
const sans = "'Outfit', sans-serif";

// ============================================================
// SMALL COMPONENTS
// ============================================================

function Field({label, value, onChange, min, max, step=0.1, unit, width}) {
  return (
    <div style={{marginBottom:8, width: width||'100%'}}>
      <label style={{display:'block',fontSize:10,color:C.dim,marginBottom:2}}>
        {label} {unit && <span style={{color:C.bd2}}>({unit})</span>}
      </label>
      <input type="number" value={value} min={min} max={max} step={step}
        onChange={e=>onChange(parseFloat(e.target.value)||0)}
        style={{
          width:'100%',padding:'5px 8px',borderRadius:5,background:C.bg,
          border:`1px solid ${C.bd}`,color:C.tx,fontFamily:mono,fontSize:12,outline:'none'
        }}
      />
    </div>
  );
}

function Card({label,value,unit,color,sub}) {
  return (
    <div style={{background:C.sf2,borderRadius:8,padding:'10px 14px',border:`1px solid ${C.bd}`,flex:'1 1 120px',minWidth:115}}>
      <div style={{fontSize:10,color:C.dim}}>{label}</div>
      <div style={{fontSize:20,fontWeight:700,color:color||C.tx,fontFamily:mono}}>
        {typeof value==='number'?value.toFixed(2):value}
        <span style={{fontSize:11,fontWeight:400,marginLeft:3,color:C.dim}}>{unit}</span>
      </div>
      {sub&&<div style={{fontSize:9,color:C.dim,marginTop:1}}>{sub}</div>}
    </div>
  );
}

function MatDisp({m, label, unit}) {
  return (
    <div style={{marginBottom:14}}>
      <div style={{color:C.ac,fontWeight:700,marginBottom:4,fontFamily:mono,fontSize:12}}>
        [{label}] <span style={{color:C.dim,fontWeight:400}}>({unit})</span>
      </div>
      <div style={{display:'grid',gridTemplateColumns:'repeat(3,1fr)',gap:2,background:C.bg,borderRadius:6,padding:5,border:`1px solid ${C.bd}`}}>
        {m.flat().map((v,i)=>(
          <div key={i} style={{textAlign:'right',padding:'3px 6px',fontFamily:mono,fontSize:10,
            color:Math.abs(v)<0.01?C.dim:C.tx,borderRadius:3,
            background:i%4===0?'rgba(56,189,248,0.04)':'transparent'
          }}>{Math.abs(v)<0.001?'≈0':v.toFixed(2)}</div>
        ))}
      </div>
    </div>
  );
}

function LayupVis({angles, t}) {
  const n=angles.length, h=n*t;
  const bh=Math.max(6,Math.min(18, 280/n));
  return (
    <div style={{display:'flex',flexDirection:'column',gap:1}}>
      <div style={{display:'flex',justifyContent:'space-between',fontSize:9,color:C.dim,marginBottom:3}}>
        <span>Top (z=+{(h/2).toFixed(1)}mm)</span><span>Ply | Angle</span>
      </div>
      {[...angles].reverse().map((a,i)=>{
        const idx=n-1-i;
        return (
          <div key={i} style={{display:'flex',alignItems:'center',gap:6,height:bh}}>
            <div style={{flex:1,height:'100%',borderRadius:2,background:plyColor(a),opacity:0.75,position:'relative',overflow:'hidden'}}>
              <svg width="100%" height="100%" style={{position:'absolute'}}>
                {Array.from({length:10},(_,l)=>{
                  const sp=100/10,rad=a*Math.PI/180,dx=Math.cos(rad)*40,dy=Math.sin(rad)*40;
                  return <line key={l} x1={`${l*sp-dx}%`} y1={`${50-dy}%`} x2={`${l*sp+dx}%`} y2={`${50+dy}%`} stroke="rgba(0,0,0,0.2)" strokeWidth="1"/>;
                })}
              </svg>
            </div>
            <div style={{fontFamily:mono,fontSize:9,color:C.tx,minWidth:55,textAlign:'right'}}>
              {idx+1}|{a>=0?'+':''}{a}°
            </div>
          </div>
        );
      })}
      <div style={{display:'flex',justifyContent:'space-between',fontSize:9,color:C.dim,marginTop:3}}>
        <span>Bot (z=−{(h/2).toFixed(1)}mm)</span><span style={{color:C.ac}}>↕ mid z=0</span>
      </div>
    </div>
  );
}

function StressPlot({data, angles, comp, label, unit}) {
  if(!data||!data.length) return null;
  let mn=Infinity,mx=-Infinity;
  const pts=[];
  data.forEach((p)=>{
    p.pts.forEach(pt=>{
      const v=pt.stress[comp]; pts.push({z:pt.z,v});
      if(v<mn)mn=v; if(v>mx)mx=v;
    });
  });
  const pad=Math.max(Math.abs(mx),Math.abs(mn))*0.15||1;
  const xMin=Math.min(mn-pad,-pad*0.3), xMax=Math.max(mx+pad,pad*0.3);
  const zMin=pts[0]?.z??-2, zMax=pts[pts.length-1]?.z??2;
  const W=360,H=260,p={t:25,r:25,b:35,l:50};
  const pw=W-p.l-p.r, ph=H-p.t-p.b;
  const sx=v=>p.l+((v-xMin)/(xMax-xMin))*pw;
  const sy=z=>p.t+ph-((z-zMin)/(zMax-zMin))*ph;
  let pathD='';
  for(let k=0;k<data.length;k++){
    const b=data[k].pts[0],t2=data[k].pts[2];
    if(k===0) pathD+=`M ${sx(b.stress[comp])} ${sy(b.z)}`;
    else pathD+=` L ${sx(b.stress[comp])} ${sy(b.z)}`;
    pathD+=` L ${sx(t2.stress[comp])} ${sy(t2.z)}`;
  }
  const colors=['#38bdf8','#f472b6','#34d399'];
  const col=colors[comp]||'#8896b0';
  return (
    <div>
      <div style={{color:C.ac2,fontWeight:600,fontSize:12,marginBottom:3}}>{label} ({unit})</div>
      <svg width={W} height={H} style={{background:C.bg,borderRadius:7,border:`1px solid ${C.bd}`}}>
        <line x1={sx(0)} y1={p.t} x2={sx(0)} y2={p.t+ph} stroke={C.bd2} strokeWidth="1" strokeDasharray="4 2"/>
        <line x1={p.l} y1={sy(0)} x2={p.l+pw} y2={sy(0)} stroke={C.ac} strokeWidth="0.8" strokeDasharray="3 2" opacity="0.4"/>
        <text x={p.l+pw+2} y={sy(0)+3} fill={C.ac} fontSize="8" opacity="0.6">mid</text>
        {data.map((pl,k)=><line key={k} x1={p.l} y1={sy(pl.zB)} x2={p.l+pw} y2={sy(pl.zB)} stroke={C.bd} strokeWidth="0.4" opacity="0.3"/>)}
        <path d={pathD+` L ${sx(0)} ${sy(zMax)} L ${sx(0)} ${sy(zMin)} Z`} fill={col} opacity="0.1"/>
        <path d={pathD} fill="none" stroke={col} strokeWidth="2"/>
        {[0,.25,.5,.75,1].map(f=>{const v=xMin+f*(xMax-xMin);
          return <text key={f} x={sx(v)} y={p.t+ph+13} fill={C.dim} fontSize="8" textAnchor="middle">{v.toFixed(0)}</text>;})}
        <text x={p.l+pw/2} y={H-3} fill={C.dim} fontSize="9" textAnchor="middle">{label} ({unit})</text>
        {[zMin,zMin/2,0,zMax/2,zMax].map((zv,i)=>
          <text key={i} x={p.l-5} y={sy(zv)+3} fill={C.dim} fontSize="8" textAnchor="end">{zv.toFixed(1)}</text>)}
        <text x={10} y={p.t+ph/2} fill={C.dim} fontSize="9" textAnchor="middle" transform={`rotate(-90,10,${p.t+ph/2})`}>z (mm)</text>
      </svg>
    </div>
  );
}

// ============================================================
// MAIN APP
// ============================================================

const DEF_FIBER = { Ef1:231, Ef2:15, Gf12:15, nuf12:0.2, name:'AS4 Carbon Fiber' };
const DEF_MATRIX = { Em:4.1, Gm:1.5, num:0.38, name:'PEKK (neat)' };
const DEF_MATRIX_SCF = { Em:8.0, Gm:3.0, num:0.35, name:'SCF-PEKK (8.1wt% SCF)' };

export default function App() {
  const [fiber, setFiber] = useState(DEF_FIBER);
  const [matrix, setMatrix] = useState(DEF_MATRIX_SCF);
  const [Vf, setVf] = useState(0.55);
  const [plyT, setPlyT] = useState(0.2);
  const [layupText, setLayupText] = useState('0/90/45/-45');
  const [repeats, setRepeats] = useState(5);
  const [Mx, setMx] = useState(1);
  const [My, setMy] = useState(0);
  const [Mxy, setMxy] = useState(0);
  const [Nx, setNx] = useState(0);
  const [Ny, setNy] = useState(0);
  const [Nxy, setNxy] = useState(0);
  const [tab, setTab] = useState('overview');
  const [matPreset, setMatPreset] = useState('scf-pekk');

  const updF = useCallback((k,v)=>setFiber(p=>({...p,[k]:v})),[]);
  const updM = useCallback((k,v)=>setMatrix(p=>({...p,[k]:v})),[]);

  const applyPreset = (preset) => {
    setMatPreset(preset);
    if(preset==='pekk-neat') setMatrix(DEF_MATRIX);
    else if(preset==='scf-pekk') setMatrix(DEF_MATRIX_SCF);
    else if(preset==='epoxy') setMatrix({Em:3.5, Gm:1.3, num:0.35, name:'Epoxy (generic)'});
    else if(preset==='peek') setMatrix({Em:4.0, Gm:1.4, num:0.37, name:'PEEK (neat)'});
  };

  const layup = useMemo(()=>{
    const base = layupText.split('/').map(s=>parseFloat(s.trim())).filter(n=>!isNaN(n));
    const full=[]; for(let r=0;r<repeats;r++) full.push(...base);
    return full;
  },[layupText,repeats]);

  const ply = useMemo(()=>micromechanics(fiber,matrix,Vf),[fiber,matrix,Vf]);

  const results = useMemo(()=>{
    if(!layup.length) return null;
    const Q = computeQ(ply.E1, ply.E2, ply.G12, ply.nu12);
    const abd = computeABD(layup, plyT, Q);
    const ip = effInplane(abd.A, abd.h);
    const fl = effFlexural(abd.D, abd.h);
    const {eps0,kappa} = solveResponse(abd.A,abd.B,abd.D,[Nx,Ny,Nxy],[Mx,My,Mxy]);
    const stresses = plyStresses(abd.qbars, abd.z, eps0, kappa, layup);
    return {Q,abd,ip,fl,eps0,kappa,stresses};
  },[layup,plyT,ply,Nx,Ny,Nxy,Mx,My,Mxy]);

  const tabs = [
    {id:'overview',label:'Overview'},
    {id:'micro',label:'Micromechanics'},
    {id:'abd',label:'ABD Matrices'},
    {id:'stress',label:'Stress Profile'},
    {id:'glossary',label:'Theory Guide'},
  ];

  const sectionHead = (text,color) => (
    <div style={{fontSize:12,fontWeight:700,color:color||C.ac,marginBottom:10,marginTop:14,
      paddingBottom:4,borderBottom:`1px solid ${C.bd}`,letterSpacing:'0.3px'}}>{text}</div>
  );

  return (
    <div style={{minHeight:'100vh',background:C.bg,color:C.tx,fontFamily:sans}}>
      <link href="https://fonts.googleapis.com/css2?family=Outfit:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500;700&display=swap" rel="stylesheet"/>
      {/* Header */}
      <div style={{padding:'16px 20px',borderBottom:`1px solid ${C.bd}`,
        background:'linear-gradient(135deg,rgba(56,189,248,0.06),rgba(129,140,248,0.04))'}}>
        <h1 style={{margin:0,fontSize:20,fontWeight:700,letterSpacing:'-0.5px'}}>
          <span style={{color:C.ac}}>CLT</span> Composite Laminate Analyzer
        </h1>
        <p style={{margin:'3px 0 0',fontSize:11,color:C.dim}}>
          Fiber + Matrix → Micromechanics → Classical Lamination Theory
        </p>
      </div>

      <div style={{display:'flex',flexWrap:'wrap'}}>
        {/* ── LEFT PANEL: INPUTS ── */}
        <div style={{width:290,padding:'14px 16px',borderRight:`1px solid ${C.bd}`,background:C.sf,
          flexShrink:0,maxHeight:'calc(100vh - 70px)',overflowY:'auto'}}>

          {/* FIBER */}
          {sectionHead('🧵 Fiber Properties','#38bdf8')}
          <div style={{fontSize:10,color:C.dim,marginBottom:6,fontStyle:'italic'}}>Pre-filled: AS4 (from your datasheet)</div>
          <Field label="Ef₁ — Longitudinal Modulus" unit="GPa" value={fiber.Ef1} onChange={v=>updF('Ef1',v)} step={1}/>
          <Field label="Ef₂ — Transverse Modulus" unit="GPa" value={fiber.Ef2} onChange={v=>updF('Ef2',v)} step={0.5}/>
          <Field label="Gf₁₂ — Shear Modulus" unit="GPa" value={fiber.Gf12} onChange={v=>updF('Gf12',v)} step={0.5}/>
          <Field label="νf₁₂ — Poisson's Ratio" value={fiber.nuf12} onChange={v=>updF('nuf12',v)} step={0.01}/>

          {/* MATRIX */}
          {sectionHead('🧪 Matrix Properties','#818cf8')}
          <div style={{display:'flex',gap:4,marginBottom:8,flexWrap:'wrap'}}>
            {[['scf-pekk','SCF-PEKK'],['pekk-neat','PEKK neat'],['peek','PEEK'],['epoxy','Epoxy']].map(([id,label])=>(
              <button key={id} onClick={()=>applyPreset(id)} style={{
                padding:'3px 8px',borderRadius:4,border:`1px solid ${matPreset===id?C.ac:C.bd}`,
                background:matPreset===id?'rgba(56,189,248,0.15)':'transparent',
                color:matPreset===id?C.ac:C.dim,fontSize:10,cursor:'pointer',fontFamily:'inherit'
              }}>{label}</button>
            ))}
          </div>
          <Field label="Em — Matrix Modulus" unit="GPa" value={matrix.Em} onChange={v=>updM('Em',v)} step={0.1}/>
          <Field label="Gm — Matrix Shear Modulus" unit="GPa" value={matrix.Gm} onChange={v=>updM('Gm',v)} step={0.1}/>
          <Field label="νm — Matrix Poisson's Ratio" value={matrix.num} onChange={v=>updM('num',v)} step={0.01}/>

          {/* Vf */}
          {sectionHead('📐 Volume Fraction & Layup','#34d399')}
          <div style={{marginBottom:8}}>
            <label style={{display:'block',fontSize:10,color:C.dim,marginBottom:2}}>
              Vf — Fiber Volume Fraction: <b style={{color:C.gn}}>{(Vf*100).toFixed(0)}%</b>
            </label>
            <input type="range" min="0.2" max="0.75" step="0.01" value={Vf} onChange={e=>setVf(parseFloat(e.target.value))}
              style={{width:'100%',accentColor:C.gn}}/>
            <div style={{display:'flex',justifyContent:'space-between',fontSize:9,color:C.dim}}>
              <span>20%</span><span>75%</span>
            </div>
          </div>
          <Field label="Ply Thickness" unit="mm" value={plyT} onChange={setPlyT} step={0.01}/>
          <div style={{marginBottom:8}}>
            <label style={{display:'block',fontSize:10,color:C.dim,marginBottom:2}}>Base Layup (e.g. 0/90/45/-45)</label>
            <input value={layupText} onChange={e=>setLayupText(e.target.value)}
              style={{width:'100%',padding:'5px 8px',borderRadius:5,background:C.bg,
                border:`1px solid ${C.bd}`,color:C.tx,fontFamily:mono,fontSize:12,outline:'none'}}/>
          </div>
          <Field label="Repeats" value={repeats} onChange={v=>setRepeats(Math.max(1,Math.round(v)))} step={1}/>
          <div style={{fontSize:10,color:C.dim}}>
            Total: {layup.length} plies = <b style={{color:C.tx}}>{(layup.length*plyT).toFixed(2)} mm</b>
          </div>

          {/* Loading */}
          {sectionHead('⚡ Applied Loading','#fb923c')}
          <div style={{fontSize:9,color:C.dim,marginBottom:6}}>Force (N/mm) & Moment (N·mm/mm)</div>
          <div style={{display:'grid',gridTemplateColumns:'1fr 1fr 1fr',gap:4}}>
            <Field label="Nₓ" value={Nx} onChange={setNx} step={1}/>
            <Field label="Nᵧ" value={Ny} onChange={setNy} step={1}/>
            <Field label="Nₓᵧ" value={Nxy} onChange={setNxy} step={1}/>
            <Field label="Mₓ" value={Mx} onChange={setMx} step={0.1}/>
            <Field label="Mᵧ" value={My} onChange={setMy} step={0.1}/>
            <Field label="Mₓᵧ" value={Mxy} onChange={setMxy} step={0.1}/>
          </div>
        </div>

        {/* ── RIGHT PANEL: RESULTS ── */}
        <div style={{flex:1,minWidth:0}}>
          <div style={{display:'flex',borderBottom:`1px solid ${C.bd}`,background:C.sf,overflowX:'auto'}}>
            {tabs.map(t=>(
              <button key={t.id} onClick={()=>setTab(t.id)} style={{
                padding:'9px 16px',border:'none',cursor:'pointer',whiteSpace:'nowrap',
                background:tab===t.id?C.sf2:'transparent',
                color:tab===t.id?C.ac:C.dim,fontWeight:tab===t.id?600:400,fontSize:12,
                borderBottom:tab===t.id?`2px solid ${C.ac}`:'2px solid transparent',
                fontFamily:'inherit',transition:'all 0.15s'
              }}>{t.label}</button>
            ))}
          </div>

          <div style={{padding:20,maxHeight:'calc(100vh - 115px)',overflowY:'auto'}}>

            {/* ── OVERVIEW TAB ── */}
            {tab==='overview' && results && (
              <div>
                {/* Computed ply props banner */}
                <div style={{background:'rgba(52,211,153,0.06)',borderRadius:8,padding:14,
                  border:`1px solid rgba(52,211,153,0.2)`,marginBottom:18}}>
                  <div style={{fontSize:12,fontWeight:600,color:C.gn,marginBottom:6}}>
                    Computed Ply Properties (from Micromechanics, Vf={( Vf*100).toFixed(0)}%)
                  </div>
                  <div style={{display:'flex',flexWrap:'wrap',gap:8}}>
                    <Card label="E₁" value={ply.E1} unit="GPa" color={C.ac}/>
                    <Card label="E₂" value={ply.E2} unit="GPa" color={C.ac} sub={`ROM: ${ply.E2_rom.toFixed(1)}`}/>
                    <Card label="G₁₂" value={ply.G12} unit="GPa" color={C.gn} sub={`ROM: ${ply.G12_rom.toFixed(1)}`}/>
                    <Card label="ν₁₂" value={ply.nu12} unit="" color={C.ac2}/>
                  </div>
                  <div style={{fontSize:9,color:C.dim,marginTop:6}}>
                    E₂ & G₁₂ use Halpin-Tsai model (more accurate than inverse ROM shown below each)
                  </div>
                </div>

                {/* Layup */}
                <h3 style={{fontSize:14,fontWeight:600,margin:'0 0 10px'}}>
                  Laminate: [{layupText}]<sub>{repeats}</sub> — {layup.length} plies, {(layup.length*plyT).toFixed(1)} mm
                </h3>
                <div style={{display:'flex',gap:12,marginBottom:10}}>
                  {[0,90,45,-45].map(a=>(
                    <div key={a} style={{display:'flex',alignItems:'center',gap:5,fontSize:10}}>
                      <div style={{width:10,height:10,borderRadius:2,background:plyColor(a)}}/> <span style={{color:C.dim}}>{a>=0?'+':''}{a}°</span>
                    </div>
                  ))}
                </div>
                <LayupVis angles={layup} t={plyT}/>

                {/* In-plane */}
                <h3 style={{fontSize:13,fontWeight:600,margin:'20px 0 10px',color:C.ac}}>In-Plane Properties (from [A])</h3>
                <div style={{display:'flex',flexWrap:'wrap',gap:8,marginBottom:16}}>
                  <Card label="Eₓ" value={results.ip?.Ex} unit="GPa" color={C.ac}/>
                  <Card label="Eᵧ" value={results.ip?.Ey} unit="GPa" color={C.ac}/>
                  <Card label="Gₓᵧ" value={results.ip?.Gxy} unit="GPa" color={C.gn}/>
                  <Card label="νₓᵧ" value={results.ip?.nuxy} unit="" color={C.ac2}/>
                </div>

                {/* Flexural */}
                <h3 style={{fontSize:13,fontWeight:600,margin:'20px 0 10px',color:C.pk}}>Flexural Properties (from [D]) — 4pt Bend Measures This</h3>
                <div style={{display:'flex',flexWrap:'wrap',gap:8,marginBottom:16}}>
                  <Card label="Eᶠₓ" value={results.fl?.Efx} unit="GPa" color={C.pk} sub="Paper: 53.2±3.4 GPa"/>
                  <Card label="Eᶠᵧ" value={results.fl?.Efy} unit="GPa" color={C.pk}/>
                  <Card label="Gᶠₓᵧ" value={results.fl?.Gfxy} unit="GPa" color={C.gn}/>
                  <Card label="νᶠₓᵧ" value={results.fl?.nufxy} unit="" color={C.ac2}/>
                </div>

                {/* Paper comparison */}
                <div style={{background:'rgba(56,189,248,0.06)',borderRadius:8,padding:14,border:`1px solid rgba(56,189,248,0.2)`,marginBottom:16}}>
                  <div style={{fontSize:12,fontWeight:600,color:C.ac,marginBottom:6}}>Comparison with Sharma et al. (2025)</div>
                  <table style={{width:'100%',fontSize:11,color:C.dim,borderCollapse:'collapse'}}>
                    <thead>
                      <tr style={{borderBottom:`1px solid ${C.bd}`}}>
                        <th style={{textAlign:'left',padding:'4px 8px',fontWeight:600,color:C.tx}}>Layup</th>
                        <th style={{textAlign:'right',padding:'4px 8px',fontWeight:600,color:C.tx}}>Paper Eᶠ</th>
                        <th style={{textAlign:'right',padding:'4px 8px',fontWeight:600,color:C.tx}}>CLT Eᶠ</th>
                      </tr>
                    </thead>
                    <tbody>
                      {[['[0/90]',67.6,'±10.3'],['[±45]',16.2,'±1.1'],['[0/90/±45]',53.2,'±3.4']].map(([n,v,sd],i)=>{
                        let cltV = '—';
                        if(i===0){const a2=computeABD([0,90,0,90,0,90,0,90,0,90,0,90,0,90,0,90,0,90,0,90],plyT,results.Q);
                          cltV=effFlexural(a2.D,a2.h)?.Efx?.toFixed(1);}
                        else if(i===1){const a2=computeABD([45,-45,45,-45,45,-45,45,-45,45,-45,45,-45,45,-45,45,-45,45,-45,45,-45],plyT,results.Q);
                          cltV=effFlexural(a2.D,a2.h)?.Efx?.toFixed(1);}
                        else { cltV=results.fl?.Efx?.toFixed(1); }
                        return (
                          <tr key={i} style={{borderBottom:`1px solid ${C.bd}`}}>
                            <td style={{padding:'4px 8px',fontFamily:mono,color:C.tx}}>{n}</td>
                            <td style={{padding:'4px 8px',textAlign:'right'}}>{v} {sd} GPa</td>
                            <td style={{padding:'4px 8px',textAlign:'right',color:C.ac,fontWeight:600}}>{cltV} GPa</td>
                          </tr>
                        );
                      })}
                    </tbody>
                  </table>
                </div>

                {/* B coupling */}
                <div style={{background:C.sf2,borderRadius:8,padding:14,border:`1px solid ${C.bd}`}}>
                  <div style={{fontSize:12,fontWeight:600,color:C.or,marginBottom:6}}>Coupling [B] Check</div>
                  <div style={{fontSize:11,color:C.dim}}>
                    {results.abd.B.flat().every(v=>Math.abs(v)<0.01) ?
                      <span style={{color:C.gn}}>✓ [B] ≈ 0 — Symmetric, no coupling</span> :
                      <span style={{color:C.or}}>⚠ [B] ≠ 0 — Bending-stretching coupling (non-symmetric layup). Consider [0/90/±45]₂ₛ for symmetric.</span>}
                  </div>
                </div>
              </div>
            )}

            {/* ── MICROMECHANICS TAB ── */}
            {tab==='micro' && (
              <div style={{maxWidth:700}}>
                <h3 style={{fontSize:15,fontWeight:600,margin:'0 0 14px'}}>Micromechanics: Fiber + Matrix → Ply</h3>

                <div style={{background:C.sf2,borderRadius:8,padding:16,border:`1px solid ${C.bd}`,marginBottom:18}}>
                  <div style={{fontSize:12,fontWeight:600,color:C.ac,marginBottom:10}}>Input Summary</div>
                  <div style={{display:'grid',gridTemplateColumns:'1fr 1fr',gap:16,fontSize:11}}>
                    <div>
                      <div style={{fontWeight:600,color:C.gn,marginBottom:4}}>Fiber: {fiber.name||'Custom'}</div>
                      <div style={{color:C.dim,lineHeight:1.8,fontFamily:mono,fontSize:10}}>
                        Ef₁ = {fiber.Ef1} GPa<br/>Ef₂ = {fiber.Ef2} GPa<br/>
                        Gf₁₂ = {fiber.Gf12} GPa<br/>νf₁₂ = {fiber.nuf12}
                      </div>
                    </div>
                    <div>
                      <div style={{fontWeight:600,color:C.ac2,marginBottom:4}}>Matrix: {matrix.name||'Custom'}</div>
                      <div style={{color:C.dim,lineHeight:1.8,fontFamily:mono,fontSize:10}}>
                        Em = {matrix.Em} GPa<br/>Gm = {matrix.Gm} GPa<br/>νm = {matrix.num}
                      </div>
                    </div>
                  </div>
                  <div style={{marginTop:10,fontSize:11,color:C.yl,fontWeight:600}}>Vf = {(Vf*100).toFixed(0)}%, Vm = {((1-Vf)*100).toFixed(0)}%</div>
                </div>

                <div style={{fontSize:12,fontWeight:600,color:C.gn,marginBottom:10}}>Computed Ply Properties</div>
                <table style={{width:'100%',borderCollapse:'collapse',marginBottom:20}}>
                  <thead>
                    <tr style={{borderBottom:`2px solid ${C.bd}`}}>
                      <th style={{textAlign:'left',padding:6,color:C.tx,fontSize:11}}>Property</th>
                      <th style={{textAlign:'right',padding:6,color:C.tx,fontSize:11}}>Halpin-Tsai</th>
                      <th style={{textAlign:'right',padding:6,color:C.tx,fontSize:11}}>Inv. ROM</th>
                      <th style={{textAlign:'left',padding:6,color:C.tx,fontSize:11}}>Method</th>
                    </tr>
                  </thead>
                  <tbody>
                    {[
                      ['E₁',ply.E1.toFixed(2),ply.E1.toFixed(2),'Rule of Mixtures (Voigt)'],
                      ['E₂',ply.E2.toFixed(2),ply.E2_rom.toFixed(2),'Halpin-Tsai (ξ=2)'],
                      ['G₁₂',ply.G12.toFixed(2),ply.G12_rom.toFixed(2),'Halpin-Tsai (ξ=1)'],
                      ['ν₁₂',ply.nu12.toFixed(4),'—','Rule of Mixtures'],
                      ['ν₂₁',ply.nu21.toFixed(4),'—','Reciprocal: ν₁₂·E₂/E₁'],
                    ].map(([p,ht,rom,method],i)=>(
                      <tr key={i} style={{borderBottom:`1px solid ${C.bd}`}}>
                        <td style={{padding:6,fontFamily:mono,fontSize:11,color:C.ac}}>{p}</td>
                        <td style={{padding:6,textAlign:'right',fontFamily:mono,fontSize:11,fontWeight:600,color:C.tx}}>{ht} GPa</td>
                        <td style={{padding:6,textAlign:'right',fontFamily:mono,fontSize:11,color:C.dim}}>{rom} GPa</td>
                        <td style={{padding:6,fontSize:10,color:C.dim}}>{method}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>

                {/* Equations */}
                <div style={{background:C.sf2,borderRadius:8,padding:16,border:`1px solid ${C.bd}`,fontSize:11,fontFamily:mono,lineHeight:2.2,color:C.dim}}>
                  <div style={{fontWeight:600,color:C.tx,marginBottom:6,fontFamily:sans}}>Micromechanics Equations Used</div>
                  <span style={{color:C.ac}}>E₁</span> = Vf·Ef₁ + Vm·Em = {Vf.toFixed(2)}×{fiber.Ef1} + {(1-Vf).toFixed(2)}×{matrix.Em} = <b style={{color:C.tx}}>{ply.E1.toFixed(2)} GPa</b><br/>
                  <span style={{color:C.ac}}>E₂</span> = Em·(1+ξη·Vf)/(1−η·Vf), η=(Ef₂/Em−1)/(Ef₂/Em+ξ), ξ=2 → <b style={{color:C.tx}}>{ply.E2.toFixed(2)} GPa</b><br/>
                  <span style={{color:C.ac}}>G₁₂</span> = Gm·(1+ξη·Vf)/(1−η·Vf), η=(Gf/Gm−1)/(Gf/Gm+ξ), ξ=1 → <b style={{color:C.tx}}>{ply.G12.toFixed(2)} GPa</b><br/>
                  <span style={{color:C.ac}}>ν₁₂</span> = Vf·νf + Vm·νm = {Vf.toFixed(2)}×{fiber.nuf12} + {(1-Vf).toFixed(2)}×{matrix.num} = <b style={{color:C.tx}}>{ply.nu12.toFixed(4)}</b>
                </div>
              </div>
            )}

            {/* ── ABD TAB ── */}
            {tab==='abd' && results && (
              <div>
                <h3 style={{fontSize:14,fontWeight:600,margin:'0 0 14px'}}>ABD Stiffness Matrices</h3>
                <div style={{display:'flex',flexWrap:'wrap',gap:16}}>
                  <MatDisp m={results.abd.A} label="A" unit="GPa·mm"/>
                  <MatDisp m={results.abd.B} label="B" unit="GPa·mm²"/>
                  <MatDisp m={results.abd.D} label="D" unit="GPa·mm³"/>
                </div>
                <div style={{marginTop:20}}>
                  <MatDisp m={[[results.Q.Q11,results.Q.Q12,0],[results.Q.Q12,results.Q.Q22,0],[0,0,results.Q.Q66]]} label="Q" unit="GPa — ply level"/>
                  <div style={{fontSize:10,color:C.dim}}>ν₂₁ = {results.Q.nu21.toFixed(5)}</div>
                </div>
              </div>
            )}

            {/* ── STRESS TAB ── */}
            {tab==='stress' && results && (
              <div>
                <h3 style={{fontSize:14,fontWeight:600,margin:'0 0 6px'}}>Through-Thickness Stress Distribution</h3>
                <div style={{fontSize:10,color:C.dim,marginBottom:14}}>
                  N=[{Nx},{Ny},{Nxy}] N/mm | M=[{Mx},{My},{Mxy}] N·mm/mm
                </div>
                <div style={{display:'flex',flexWrap:'wrap',gap:6,marginBottom:16}}>
                  <Card label="εₓ⁰" value={(results.eps0[0]*1e6).toFixed(1)} unit="µε" color={C.ac}/>
                  <Card label="εᵧ⁰" value={(results.eps0[1]*1e6).toFixed(1)} unit="µε" color={C.pk}/>
                  <Card label="κₓ" value={results.kappa[0].toExponential(3)} unit="1/mm" color={C.or}/>
                </div>
                <div style={{display:'flex',flexWrap:'wrap',gap:16}}>
                  <StressPlot data={results.stresses} angles={layup} comp={0} label="σₓ" unit="MPa"/>
                  <StressPlot data={results.stresses} angles={layup} comp={1} label="σᵧ" unit="MPa"/>
                  <StressPlot data={results.stresses} angles={layup} comp={2} label="τₓᵧ" unit="MPa"/>
                </div>
                <div style={{background:C.sf2,borderRadius:8,padding:14,marginTop:16,border:`1px solid ${C.bd}`,fontSize:11,color:C.dim,lineHeight:1.8}}>
                  <b style={{color:C.tx}}>Reading stress plots:</b><br/>
                  • Vertical axis = z position (bottom to top of laminate)<br/>
                  • Under bending (Mₓ), top plies compress, bottom plies stretch (or vice versa)<br/>
                  • Jumps at ply boundaries = each ply angle has different Q̄ in global frame<br/>
                  • ±45° plies show different σₓ than 0°/90° plies because their transformed stiffness differs
                </div>

                {/* Ply-by-ply table */}
                <h3 style={{fontSize:13,fontWeight:600,margin:'20px 0 8px',color:C.ac2}}>Ply-by-Ply Stresses (at mid-ply)</h3>
                <div style={{overflowX:'auto'}}>
                  <table style={{width:'100%',borderCollapse:'collapse',fontSize:10,fontFamily:mono}}>
                    <thead>
                      <tr style={{borderBottom:`2px solid ${C.bd}`}}>
                        {['Ply','Angle','z_mid','σₓ','σᵧ','τₓᵧ'].map(h=>(
                          <th key={h} style={{padding:'4px 6px',textAlign:h==='Ply'||h==='Angle'?'center':'right',color:C.tx,fontSize:9}}>{h}</th>
                        ))}
                      </tr>
                    </thead>
                    <tbody>
                      {results.stresses.map((p,i)=>{
                        const s=p.pts[1].stress;
                        const zm=p.pts[1].z;
                        return (
                          <tr key={i} style={{borderBottom:`1px solid ${C.bd}`,background:i%2===0?'transparent':'rgba(255,255,255,0.01)'}}>
                            <td style={{padding:'3px 6px',textAlign:'center',color:C.dim}}>{i+1}</td>
                            <td style={{padding:'3px 6px',textAlign:'center',color:plyColor(p.angle)}}>{p.angle>=0?'+':''}{p.angle}°</td>
                            <td style={{padding:'3px 6px',textAlign:'right',color:C.dim}}>{zm.toFixed(2)}</td>
                            <td style={{padding:'3px 6px',textAlign:'right',color:s[0]>=0?'#38bdf8':'#f87171'}}>{s[0].toFixed(1)}</td>
                            <td style={{padding:'3px 6px',textAlign:'right',color:s[1]>=0?'#38bdf8':'#f87171'}}>{s[1].toFixed(1)}</td>
                            <td style={{padding:'3px 6px',textAlign:'right',color:C.tx}}>{s[2].toFixed(1)}</td>
                          </tr>
                        );
                      })}
                    </tbody>
                  </table>
                </div>
              </div>
            )}

            {/* ── GLOSSARY TAB ── */}
            {tab==='glossary' && (
              <div style={{fontSize:12,lineHeight:1.8,color:C.dim,maxWidth:720}}>
                <h3 style={{color:C.ac,fontWeight:700,fontSize:16}}>Composite Theory Reference</h3>
                <table style={{width:'100%',borderCollapse:'collapse',marginTop:12}}>
                  <tbody>
                    {[
                      ['Ply / Lamina','Single layer of fibers in matrix'],
                      ['Laminate','Stack of bonded plies'],
                      ['Layup','Order & angles of plies, e.g. [0/90/±45]₅'],
                      ['Vf (Fiber Volume Fraction)','% of volume that is fiber (~50-65% typical)'],
                      ['Rule of Mixtures (Voigt)','E₁ = Vf·Ef + Vm·Em — upper bound, exact for longitudinal'],
                      ['Halpin-Tsai','Semi-empirical model for E₂ and G₁₂ — more accurate than inverse ROM'],
                      ['Inverse ROM (Reuss)','1/E₂ = Vf/Ef₂ + Vm/Em — lower bound, underestimates'],
                      ['E₁ (Longitudinal)','Stiffness along fibers — fiber dominated'],
                      ['E₂ (Transverse)','Stiffness across fibers — matrix dominated'],
                      ['G₁₂ (Shear)','Resistance to in-plane shearing'],
                      ['ν₁₂ (Poisson)','Lateral contraction ratio when loaded along fibers'],
                      ['[Q] Reduced Stiffness','3×3 stress-strain relation in ply (1-2) coordinates'],
                      ['[Q̄] Transformed','[Q] rotated by angle θ into global (x-y) coordinates'],
                      ['[A] Extensional','In-plane force ↔ strain. Units: GPa·mm'],
                      ['[B] Coupling','Bending-stretching coupling. Zero for symmetric laminates'],
                      ['[D] Bending','Moment ↔ curvature. Units: GPa·mm³. Determines flexural modulus'],
                      ['Flexural Modulus','What 4-point bending measures — from [D], not [A]!'],
                      ['Quasi-Isotropic','Equal in-plane stiffness all directions (e.g. [0/±45/90])'],
                      ['Balanced','For every +θ there is a −θ ply (A₁₆=A₂₆=0)'],
                      ['Symmetric','Mirrored about midplane → [B]=0'],
                      ['Tsai-Wu','Quadratic failure criterion in material coords'],
                      ['Stress Shielding','Implant too stiff → carries load instead of bone'],
                      ['Delamination','Ply separation — primary composite failure mode'],
                      ['FFF','Fused Filament Fabrication — your 3D printing process'],
                      ['Consolidation','Post-print compression molding to remove voids'],
                    ].map(([t,d],i)=>(
                      <tr key={i} style={{borderBottom:`1px solid ${C.bd}`}}>
                        <td style={{padding:'6px 10px 6px 0',color:C.ac,fontWeight:500,verticalAlign:'top',whiteSpace:'nowrap',fontSize:11}}>{t}</td>
                        <td style={{padding:'6px 0',fontSize:11}}>{d}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>

                <div style={{marginTop:20}}>
                  <h4 style={{color:C.tx,fontWeight:600}}>Full CLT Workflow</h4>
                  <div style={{background:C.sf2,borderRadius:8,padding:14,border:`1px solid ${C.bd}`,fontSize:11,fontFamily:mono,lineHeight:2.2}}>
                    1. Input fiber props (Ef₁, Ef₂, Gf₁₂, νf) + matrix (Em, Gm, νm)<br/>
                    2. Micromechanics → ply E₁, E₂, G₁₂, ν₁₂<br/>
                    3. Build [Q] reduced stiffness<br/>
                    4. For each angle θ → [Q̄(θ)]<br/>
                    5. Integrate through thickness → [A], [B], [D]<br/>
                    6. Invert [A] → in-plane Ex, Ey, Gxy, νxy<br/>
                    7. Invert [D] → flexural Efx, Efy (what you measure!)<br/>
                    8. Apply loads → midplane strains + curvatures<br/>
                    9. Ply-by-ply stresses in global & material coords<br/>
                    10. Failure check (Tsai-Wu, Max Stress)
                  </div>
                </div>
              </div>
            )}
          </div>
        </div>
      </div>
    </div>
  );
}
