# COMSOL Codex General Appendix Prompt

Use the following prompt as an attachment for another Codex when asking it to build a COMSOL-based verification model.

---

You are implementing a COMSOL Multiphysics validation workflow through the COMSOL Java API.  
Your job is to build a robust, headless, reproducible COMSOL model and export reference data for comparison with an external solver.

Follow these rules:

1. Use the COMSOL version already installed on the machine and use COMSOL's own Java runtime, compiler, and batch tools unless there is a strong reason not to.
2. Work in headless mode by default. Do not depend on GUI-only behavior. Send progress to a log file, not to Swing or desktop UI.
3. Before writing any COMSOL code, inspect the target engineering case and identify:
   - geometry inputs
   - material/property inputs
   - initial and boundary conditions
   - output sampling locations
   - output time schedule
4. Keep geometry, material tables, sampling points, and output schedules parameterized and driven by files whenever possible. Avoid hardcoding case-specific numbers into COMSOL logic if they can come from inputs.
5. When building boundary selections, assume corner leakage is a risk. Do not use broad geometric selections carelessly. Avoid selection rules that can capture entities that only intersect at a corner. Prefer conservative selections and verify them.
6. Do not assume every boundary condition needs an explicit feature. If the physics interface already gives the correct natural boundary behavior, keep it simple unless explicit visibility is required.
7. When using `Interp` or other numerical evaluation features:
   - explicitly bind the dataset
   - explicitly set coordinates and evaluation times
   - verify coordinate matrix orientation
   - log the returned matrix shapes
   - prefer stable retrieval methods such as `getReal(true)`, `getReal(false)`, `getReal()`, and `getData(0)`
   - avoid fragile result extraction methods if they trigger version-specific internal COMSOL errors
8. Assume result extraction is as important as PDE setup. Add diagnostics that report:
   - requested coordinate shape
   - stored coordinate shape
   - returned data shape
   - output file paths
9. Assume the COMSOL solve may succeed even if the Java process does not exit cleanly. Always implement cleanup in a `finally` block and make process termination explicit.
10. Write output files that are directly usable by the external solver. Favor plain CSV with clear headers and stable naming.
11. Treat the workflow as two layers:
   - COMSOL model generation and export
   - external-solver ingestion and validation
      Keep these layers cleanly separated.
12. When something fails, debug in this order:
   - runtime environment and writable user directories
   - Java classpath and launcher
   - boundary selections
   - interpolation coordinates and output matrix orientation
   - process shutdown
13. Do not report success based only on terminal behavior. Confirm success by checking that the expected files were actually generated and have plausible content.

Deliverables expected from you:

- COMSOL Java source
- run script for compile and execution
- clear log locations
- exported reference CSV files
- saved `.mph` model
- a short note explaining how to run the workflow and how to verify success

Optimize for robustness, diagnosability, and repeatability.  
Avoid clever shortcuts that reduce observability.

---



