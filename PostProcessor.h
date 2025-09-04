#pragma once
#include <string>

class Mesh;
class FractureNetwork;
struct FieldRegistry;

class PostProcessor {

public:
	explicit PostProcessor(std::string outdir = "results");

	void exportMatrixValue(const Mesh& mesh,
		const FieldRegistry& reg,
		double time, int step) const;

	void exportFractureValue(const Mesh& mesh,
		const FractureNetwork& frNet,
		const FieldRegistry& reg_fr,
		double time, int step) const;

	void exportMatrixValueInterpolated(const Mesh& mesh,
		const FieldRegistry& reg_prev,
		const FieldRegistry& reg_curr,
		double time, int step, double alpha) const;

	void exportFractureValueInterpolated(const Mesh& mesh,
		const FractureNetwork& frNet,
		const FieldRegistry& reg_prev_fr,
		const FieldRegistry& reg_curr_fr,
		double time, int step, double alpha) const;


private:
	std::string outdir_; //输出目录

	static std::string stepName_(int step); //辅助函数：根据时间步数生成文件名

	void appendTimeIndex_(double time, int step) const; //辅助函数：将时间和时间步数写入文件

	// 实际写文件
	void exportMatrix_(const Mesh& mesh,
		const FieldRegistry& regA,
		const FieldRegistry* regB, double alpha,
		double time, int step) const;

	void exportFracture_(const Mesh& mesh,
		const FractureNetwork& frNet,
		const FieldRegistry& regA,
		const FieldRegistry* regB, double alpha,
		double time, int step) const;

};