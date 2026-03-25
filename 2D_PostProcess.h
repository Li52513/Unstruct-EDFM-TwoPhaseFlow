#pragma once

#include <string>
#include <vector>
#include <memory>
#include <set>
#include <algorithm>

// ����������� (���� 2D ����)
#include "MeshManager.h" 
#include "2D_FieldManager.h"
#include "VTK_BoundaryVisualization.h"

/**
 * @class PostProcess_2D
 * @brief 2D-EDFM �������ӻ�ģ�� (Tecplot Exporter)
 * @details �� 2D ����(Triangle/Quad) �� 1D �ѷ�(LineSeg) ��������������Ϊ Tecplot ��ʽ��
 */
class PostProcess_2D
{
public:
    /**
     * @brief ���캯��
     * @param meshMgr 2D���������
     * @param fieldMgr 2D��������
     */
    PostProcess_2D(const MeshManager& meshMgr,
                   const FieldManager_2D& fieldMgr,
                   const VTKBoundaryVisualizationContext* bcVizCtx = nullptr);

    /**
     * @brief ����ȫ�����ݵ� Tecplot
     * @param filename ����ļ��� (�� "Result_2D.dat")
     * @param time ��ǰģ��ʱ��
     */
    void ExportTecplot(const std::string& filename, double time = 0.0) const;

    /**
     * @brief ����ȫ�����ݵ� ParaView (VTK ��ʽ)
     * @param filename ����ļ��� (�� "Result_2D.vtk")
     * @param time ��ǰģ��ʱ��
     * @details �Զ�ӳ����� (Tri/Quad) �� 1D �ѷ� (Line) �� VTK �� Unstructured Grid ��ʽ
     */
    void ExportVTK(const std::string& filename, double time = 0.0) const;

    /**
     * @brief Export to ParaView XML Unstructured Grid format (.vtu)
     * @param filename Output file path (e.g. "Result_2D_step0.vtu")
     * @param time     Current simulation time [s]
     * @details Writes ASCII VTK XML with CellData for all registered fields.
     *          Supports Triangle (type 5), Quad (type 9) matrix cells and
     *          Line (type 3) fracture segments.  Suitable for ParaView time-series
     *          animation when combined with ExportPVD().
     */
    void ExportVTU(const std::string& filename, double time = 0.0) const;

    /**
     * @brief Write a ParaView Data collection file (.pvd) for time-series animation
     * @param pvd_filename  Output .pvd file path
     * @param vtu_filenames Ordered list of .vtu file paths (relative or absolute)
     * @param times         Corresponding simulation times [s]
     * @details The .pvd file references all .vtu snapshots; open it in ParaView
     *          to animate the transient solution.
     */
    static void ExportPVD(const std::string& pvd_filename,
                          const std::vector<std::string>& vtu_filenames,
                          const std::vector<double>& times);

    /**
     * @brief [FIM ���ݽ�ά����] �� 2D �Զ�΢�ֳ� (ADVar) ͬ��Ϊ��˫���ȱ����� (double)
     * @tparam N �����Ա��������� (�絥��Ϊ2������Ϊ3)
     * @param fieldMgr 2D ������������ (�� const�����贴��/�޸��³�)
     * @param adFieldName ԭʼ���ݶȵ��Զ�΢�ֳ����� (���� "Pressure_AD")
     * @param outScalarName ��ά�����ڵ����Ĵ����������� (���� "Pressure")
     * @details
     * �Զ��������������ѷ��򣬰�ȫ���� ADVar.val ����ֵ�������� volScalarField��
     * ����ģ�� ExportTecplot ������� outScalarName ����ʵ����Խ����յ�����
     */
    template<int N>
    static void SyncADFieldToScalar(FieldManager_2D& fieldMgr, const std::string& adFieldName, const std::string& outScalarName)
    {
        // 1. ͬ�������� (Matrix)
        auto adMat = fieldMgr.getMatrixADScalar<N>(adFieldName);
        if (adMat) {
            auto scMat = fieldMgr.getOrCreateMatrixScalar(outScalarName);
            const size_t nCells = adMat->data.size();
            for (size_t i = 0; i < nCells; ++i) {
                scMat->data[i] = adMat->data[i].val; // ��ȷ���� ADVar ��Ա���� val
            }
        }

        // 2. ͬ���ѷ��� (Fracture)
        auto adFrac = fieldMgr.getFractureADScalar<N>(adFieldName);
        if (adFrac) {
            auto scFrac = fieldMgr.getOrCreateFractureScalar(outScalarName);
            const size_t nFracCells = adFrac->data.size();
            for (size_t i = 0; i < nFracCells; ++i) {
                scFrac->data[i] = adFrac->data[i].val; // ��ȷ���� ADVar ��Ա���� val
            }
        }
    }

private:
    const MeshManager& meshMgr_;
    const FieldManager_2D& fieldMgr_;
    const VTKBoundaryVisualizationContext* bcVizCtx_ = nullptr;

    // =========================================================
    // �ڲ���������
    // =========================================================

    /**
     * @brief ��ȡ������ (2D Matrix + 1D Fracture) ��ע���Ψһ����������
     */
    std::vector<std::string> GetAllUniqueFieldNames() const;

    /**
     * @brief ���ݽڵ����ж� Tecplot ��Ԫ����
     * @param numNodes ��Ԫ������ (2: FELINESEG, 3: FETRIANGLE, 4: FEQUADRILATERAL)
     */
    std::string GetTecplotElementType(int numNodes) const;
};
