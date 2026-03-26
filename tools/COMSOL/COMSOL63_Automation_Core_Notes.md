# COMSOL 6.3 Automation Core Notes

## Summary

This note keeps only the high-signal lessons from one full COMSOL Java API validation workflow.  
The goal is to avoid repeating the same integration mistakes when building future verification models.

```mermaid
flowchart LR
    A["Compile Java with COMSOL JRE/JDK"] --> B["Run COMSOL model headless"]
    B --> C["Solve PDE and export CSV/MPH"]
    C --> D["Read reference data in the engineering solver"]
    D --> E["Post-process and plot"]
```

## What Actually Happened

1. The engineering-side case was first run to generate COMSOL input tables, sampling points, and schedules.
2. A COMSOL 6.3 Java model was then compiled and executed in headless mode.
3. Several failures appeared during automation, but the PDE solve itself was usually not the problem.
4. The main errors were in runtime environment setup, Java launching, boundary selection, interpolation extraction, and process shutdown.
5. After fixing those issues, COMSOL successfully exported reference CSV files and `.mph` models.

## Core Lessons

### 1. Use COMSOL's own Java runtime

- Use the Java tools shipped with COMSOL whenever possible.
- Do not mix system Java with COMSOL Java API runs unless there is a very specific reason.
- Keep compile and runtime paths explicit.

### 2. Headless COMSOL needs a writable user profile

- COMSOL runtime writes under `%USERPROFILE%\.comsol\v63\...`.
- If these folders are missing or not writable, automation can fail before any real model work starts.
- Pre-create these folders in automation scripts if needed.

### 3. Headless progress logging should go to a file, not GUI

- `ModelUtil.showProgress(true)` can trigger GUI toolkit issues in headless runs.
- Prefer `ModelUtil.showProgress(<log_file_path>)`.

### 4. Keep Java classpath simple

- A giant hand-built classpath string is fragile.
- Prefer a plugin wildcard plus the script directory.
- If Java cannot find the main class, inspect classpath formatting first.

### 5. Be careful with COMSOL selections

- Box selections are not the same as explicit selections.
- Do not call APIs meant only for explicit selections on a Box selection.
- `Box + intersects` can accidentally include boundaries that only touch at corners.
- For left/right boundary picking, `Box + inside` is safer than `intersects` when corner leakage matters.

### 6. Natural boundaries already give zero flux

- In many PDE interfaces, unassigned boundaries remain natural no-flux boundaries.
- If the GUI must show them explicitly, add flux features.
- If not, avoid unnecessary extra features.

### 7. Interpolation result extraction is a major failure point

- Do not assume the `Interp` numerical feature returns data in the shape you expect.
- Explicitly bind the dataset, coordinates, and evaluation times.
- Prefer:
  - `getReal(true)`
  - `getReal(false)`
  - `getReal()`
  - `getData(0)`
- Avoid `computeResult()` if it triggers COMSOL internal property errors in your version.

### 8. Coordinate matrix orientation must be checked

- `Interp` may recognize solution times but return zero sampled points if coordinate layout is wrong.
- Log both:
  - requested coordinate matrix shape
  - stored coordinate matrix shape reported by COMSOL
- If necessary, automatically try both common layouts.

### 9. Successful solve does not guarantee clean process exit

- COMSOL can finish solving and write all output files while `java.exe` still stays alive.
- Always clean up in a `finally` block.
- A robust shutdown path is:
  - `ModelUtil.clear()`
  - `ModelUtil.disconnect()`
  - `System.exit(exitCode)`

### 10. Trust output files more than terminal silence

- A quiet terminal does not mean the run failed.
- Check whether these artifacts were produced:
  - reference CSV files
  - `.mph` model
  - progress log
  - mesh-check report
- In headless workflows, file timestamps are often the fastest truth source.

## Recommended Minimal Workflow

1. Run the engineering case first and export COMSOL input tables, sampling points, and schedules.
2. Compile Java with COMSOL's bundled Java tools.
3. Run COMSOL headless with file-based progress logging.
4. Export reference CSV files first; fine-mesh checks can be a separate second pass.
5. Re-run the engineering executable to ingest the COMSOL reference data.
6. Only then generate plots and validation summaries.

## Final Rule Set

- Parameterize everything.
- Keep selections conservative.
- Log shapes and file paths.
- Prefer stable API calls over clever ones.
- Treat interpolation and shutdown as first-class engineering tasks, not afterthoughts.

