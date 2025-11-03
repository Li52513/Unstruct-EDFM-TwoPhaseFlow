from pathlib import Path
text = Path("Solver_AssemblerCOO.h").read_text(encoding="gbk")
start = text.index('inline bool assembleTemperature_singlePhase_COO')
end = text.index(')\n{', start)
print(text[start:end+2])
