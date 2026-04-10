# CLT Composite Laminate Analyzer

**Classical Lamination Theory** analysis tool for carbon fiber / PEKK composites.

Computes ply properties from fiber + matrix inputs using micromechanics (Rule of Mixtures, Halpin-Tsai), then runs full CLT to give you ABD matrices, effective in-plane and flexural properties, and through-thickness stress distributions.

Pre-configured for **AS4 carbon fiber + SCF-PEKK matrix** with validation against published experimental data.

![screenshot](screenshot.png)

---

## Prerequisites

You need **Node.js** installed (version 18 or newer).

- **Windows**: Download from https://nodejs.org — click the LTS version, run the installer
- **Mac**: `brew install node` or download from https://nodejs.org
- **Linux**: `sudo apt install nodejs npm` or download from https://nodejs.org

To check if you have it: open a terminal and type `node --version`

---

## Setup & Run (3 commands)

Open a terminal / command prompt, navigate to this folder, then:

```bash
npm install        # downloads dependencies (one-time, takes ~30 seconds)
npm run dev        # starts the app — opens in your browser automatically
```

That's it. The app runs at **http://localhost:5173** and auto-opens in your browser.

Press `Ctrl+C` in the terminal to stop the server when done.

---

## How to Use

1. **Left panel** — Enter your fiber properties (AS4 pre-filled from datasheet), matrix properties (preset buttons for SCF-PEKK, neat PEKK, PEEK, epoxy), fiber volume fraction, and layup definition
2. **Overview tab** — See computed ply properties, layup visualization, in-plane and flexural moduli, comparison with published paper data
3. **Micromechanics tab** — Shows exactly how fiber + matrix → ply properties with Halpin-Tsai equations
4. **ABD Matrices tab** — Full [A], [B], [D] stiffness matrices and ply-level [Q] matrix
5. **Stress Profile tab** — Through-thickness σx, σy, τxy distributions with ply-by-ply stress table
6. **Theory Guide tab** — Composite terminology and CLT workflow reference

---

## Building for Sharing (no Node.js needed to open)

To create a static version anyone can open by double-clicking an HTML file:

```bash
npm run build
```

This creates a `dist/` folder. Zip it up and send it — the recipient just opens `dist/index.html` in any browser. No installation needed.

---

## Project Structure

```
composite-clt-app/
├── index.html          # HTML entry point
├── package.json        # dependencies & scripts
├── vite.config.js      # build tool config
├── README.md           # this file
└── src/
    ├── main.jsx        # React mount point
    └── App.jsx         # All analysis code (micromechanics + CLT + UI)
```

---

## Theory Reference

- **Micromechanics**: Rule of Mixtures (E₁, ν₁₂), Halpin-Tsai (E₂, G₁₂)
- **CLT**: Reduced stiffness [Q] → Transformed [Q̄] → ABD integration → Effective properties
- **Flexural modulus** from [D] matrix = what 4-point bending measures
- **In-plane modulus** from [A] matrix = tensile test result

Based on: Sharma et al. (2025) "Hybrid-manufactured silicon nitride coated CFR-PEKK" — J. Mech. Behav. Biomed. Mater.
