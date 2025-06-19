#pragma once
#include <vector>
#include <map>
#include "Node.h"
#include "Face.h"
#include "Cell.h"
#include "Matrix.h"  
#include <fstream>


class Mesh 
{
public:
 
    //------------------基岩网格信息----------------------//
       
    vector<Node> nodes;
    vector<Face> faces;
    vector<Cell> cells;
    vector<Cell> innerCells;
    vector<Cell> boundaryCells;
    map<int, Node> nodesMap;
    map<int, int> cellId2index; 
    int gridCount;
    Mesh();
    void BuildMesh(double lengthX, double lengthY, double lengthZ, int nx, int ny, int nz,bool usePrism, bool useQuadBase);
    void ClassifySolidMatrixCells();
    void printMeshInfo();
    void exportToTxt(const std::string& prefix) const;

	////------------------基岩物性参数----------------------//
 //   void assignRockProperties(const Matrix& matrix); //代替换输入参数
 //   Matrix matrixProperties;
 //   //---------------------------------------------------//
   
    
    
    double getCellArea(int cellID)  const;
    //记录边界面的ID
    vector<int> getBoundaryFaceIDs() const 
    {
        vector<int> ids;
        for (auto& f : faces)
        {
            if (f.isBoundary()) ids.push_back(f.id);
        }
        return ids;
    }
    void ComputeSolidMatrixMeshFaceGeometricInfor(NormalVectorCorrectionMethod method);

    int ghostStartIndex = 0;
    // （可选）Dirichlet 边界压力  //待优化
    void CreateSolidMatrixGhostCells();    // 新增：生成 Ghost Cell
    void exportGhostInfo(const std::string& prefix) const;

};


//#pragma once
//#include <vector>
//#include <map>
//#include <set>
//#include <algorithm>
//#include <iostream>
//#include <gmsh.h>
//#include "Node.h"
//#include "Face.h"
//#include "Cell.h"
//using namespace std;
//
//class Mesh {
//public:
//    // 全局存储所有节点、面和单元
//    vector<Node> nodes;
//    vector<Face> faces;
//    vector<Cell> cells;
//
//    // 内部单元与边界单元分别存储
//    vector<Cell> innerCells;
//    vector<Cell> boundaryCells;
//
//    int gridCount;   // 网格单元总数
//
//    // 为便于查找，建立以节点ID为key的映射，值类型为 Node 对象
//    map<int, Node> nodesMap;
//
//    Mesh() : gridCount(0) {} // 默认构造函数
//
//    // 调用 gmsh 生成矩形区域的三角形网格，并构造 Mesh 数据结构
//    void BuildMesh(double lengthX, double lengthY, int sectionNumX, int sectionNumY) 
//    {
//        gmsh::initialize();
//		gmsh::model::add("UnstructuredMesh-EDFM"); // 创建模型
//
//        // 网格尺寸控制
//        double lc = min(lengthX / sectionNumX, lengthY / sectionNumY); 
//        int p1 = gmsh::model::geo::addPoint(0, 0, 0, lc); //(x,y,z,length)
//        int p2 = gmsh::model::geo::addPoint(lengthX, 0, 0, lc);
//        int p3 = gmsh::model::geo::addPoint(0, lengthY, 0, lc);
//        int p4 = gmsh::model::geo::addPoint(lengthX, lengthY, 0, lc);
//
//        int l1 = gmsh::model::geo::addLine(p1, p2);
//        int l2 = gmsh::model::geo::addLine(p2, p4);
//        int l3 = gmsh::model::geo::addLine(p4, p3);
//        int l4 = gmsh::model::geo::addLine(p3, p1);
//
//        int curveLoop = gmsh::model::geo::addCurveLoop({ l1, l2, l3, l4 }); // 创建一个封闭的曲线环
//        int surface = gmsh::model::geo::addPlaneSurface({ curveLoop }); // 创建一个平面区域
//
//        gmsh::model::geo::synchronize(); // 同步几何模型
//		gmsh::model::mesh::generate(2); // 生成2D网格（默认三角形网格）2——2D 3——3D
//
//        // —— Step1 获取所有网格节点，并构造 Node 对象，同时建立节点ID与坐标的映射 ——// 
//        vector<size_t> nodeTags; // 节点标签
//        vector<double> nodeCoords, parametricCoords; // 节点坐标
//        gmsh::model::mesh::getNodes(nodeTags, nodeCoords, parametricCoords); // 获取所有节点信息
//        for (size_t i = 0; i < nodeTags.size(); i++) 
//        {
//            int id = static_cast<int>(nodeTags[i]);    // 将 size_t 转换为 int，得到节点编号
//            Vector coord(nodeCoords[3 * i], nodeCoords[3 * i + 1], nodeCoords[3 * i + 2]); // 读取节点坐标
//            nodes.push_back(Node(id, coord)); // 将节点信息存入 nodes 容器
//			nodesMap[id] = Node(id, coord);     // 同时建立以节点id为 key 的映射，值为 Node 对象
//        }
//
//        // —— Step2 获取三角形单元，并构造 Cell 对象 —— //
//		vector<int> elementTypes;   // 单元类型
//        vector<vector<size_t>> elementTags, nodeTagsPerElement; // 与单元类型相对应的二维数组，单元标签及单元节点标签
//        gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTagsPerElement);
//		for (size_t i = 0; i < elementTypes.size(); i++) //elementTypes的外层数组，即网格单元类型
//        {
//			if (elementTypes[i] == 2) // 2代表三角形单元 四边形单元为3
//            {  
//               
//				size_t numTriangles = elementTags[i].size(); // 三角形单元个数
//				for (size_t j = 0; j < numTriangles; j++)    // 遍历每个三角形单元
//                {
//					size_t cellId = elementTags[i][j]; // 读取单元编号
//                    vector<int> cellNodeIDs; //建立储存单元节点编号的容器	
//					// 三角形单元有三个节点
//					for (int k = 0; k < 3; k++)
//					{
//						int nodeID = static_cast<int>(nodeTagsPerElement[i][3 * j + k]); // 读取节点ID
//						cellNodeIDs.push_back(nodeID); // 将节点ID保存在cellNodeIDs容器中
//					}
//                    Cell cell(cellId, cellNodeIDs); // 构造 Cell 对象
//                    cell.computeCenterAndVolume(nodesMap); // 计算 Cell 的中心和面积
//                    cells.push_back(cell);
//                    
//                }
//            }
//        }
//        gridCount = cells.size(); // 网格单元总数
//
//        // —— 3. 构造面（Face）：利用所有单元的边构建唯一面，并建立 face 与 cell 的对应关系 —— 
//        // 使用 map 以有序边 (min, max) 作为 key 来记录共享该边的单元编号
//        map<pair<int, int>, vector<int>> faceMap;
//        for (const auto& cell : cells) 
//        {
//			const vector<int>& cn = cell.nodeIDs; // 单元的节点编号
//            if (cn.size() != 3) continue; // 仅处理三角形单元
//            // 每个三角形单元有三条边（节点编号通过 min/max 处理成有序对）
//            vector<pair<int, int>> edges = 
//            {
//                { min(cn[0], cn[1]), max(cn[0], cn[1]) },
//                { min(cn[1], cn[2]), max(cn[1], cn[2]) },
//                { min(cn[2], cn[0]), max(cn[2], cn[0]) }
//            };
//            for (auto& edge : edges) 
//            {
//				faceMap[edge].push_back(cell.id); // 将每个边节点编号排序后作为 key，将包含该边的单元编号保存在 value 中
//            }
//        }
//
//        int faceId = 1; // 面编号从1开始
//		for (auto& entry : faceMap) // 遍历 map
//        {
//            int n1 = entry.first.first;     // 访问 key 的第一个节点编号,边的第一个节点（较小的编号）
//            int n2 = entry.first.second;      // 边的第二个节点（较大的编号）
//            vector<int> faceNodeIDs = { n1, n2 };
//            vector<Vector> faceNodeCoords = { nodesMap[n1].coord, nodesMap[n2].coord };
//            Face face(faceId, faceNodeIDs, faceNodeCoords); // 构造 Face 对象
//            if (entry.second.size() == 1)
//            {
//                face.ownerCell = entry.second[0];
//                face.neighborCell = -1; // 仅被一个单元拥有，为边界面
//            }
//            else if (entry.second.size() == 2) 
//            {
//                face.ownerCell = entry.second[0];
//                face.neighborCell = entry.second[1];
//            }
//            faces.push_back(face);
//            // 同时将该面编号记录到对应单元中
//            for (int cid : entry.second) 
//            {
//                for (auto& cell : cells) 
//                {
//                    if (cell.id == cid)
//                    {
//                        cell.faceIDs.push_back(faceId);
//                        break;
//                    }
//                }
//            }
//            faceId++;
//        }
//
//        gmsh::write("UnstructuredMesh-EDFM.msh");
//        gmsh::finalize();
//
//        
//        // 目的：统一所有元素（节点、单元、面）的编号，使其均从1开始且连续
//        // 1. 重新映射节点 ID
//        map<int, int> newNodeMapping;  // 原始节点ID -> 新节点ID
//        int newId = 1;
//        for (auto& node : nodes) 
//        {
//            cout << "原始node.id为 " << node.id << endl;
//            newNodeMapping[node.id] = newId;
//            node.id = newId;
//            cout << "新node.id为 " << node.id << endl;
//            newId++;
//
//        }
//
//        // 2. 重新映射单元（Cell）的 ID
//        map<int, int> newCellMapping;  // 原始 Cell ID -> 新 Cell ID
//        newId = 1;
//        for (auto& cell : cells) 
//        {
//			cout << "原始cell.id为 " << cell.id << endl;
//            newCellMapping[cell.id] = newId;
//            cell.id = newId;
//            cout << "新cell.id为 " << cell.id << endl;
//            // 更新 cell 中的节点ID，使用新的节点映射
//            for (auto& nid : cell.nodeIDs)
//				nid = newNodeMapping[nid]; // 更新节点ID
//            newId++;
//        }
//
//        // 3. 重新映射面（Face）的 ID
//        map<int, int> newFaceMapping;  // 原始 Face ID -> 新 Face ID
//        newId = 1;
//        for (auto& face : faces)
//        {
//            cout << "原始face.id为 " << face.id << endl;
//            newFaceMapping[face.id] = newId;
//            face.id = newId;
//            cout << "新face.id为 " << face.id << endl;
//            // 更新 face 中的节点ID
//            for (auto& nid : face.nodeIDs)
//                nid = newNodeMapping[nid];
//            // 更新面引用的 Cell ID
//            if (face.ownerCell != -1)
//                face.ownerCell = newCellMapping[face.ownerCell];
//            if (face.neighborCell != -1)
//                face.neighborCell = newCellMapping[face.neighborCell];
//            newId++;
//        }
//
//        // 4. 更新每个单元中存储的面 ID（faceIDs），使用新的面映射
//        for (auto& cell : cells) 
//        {
//            for (auto& fid : cell.faceIDs)
//                fid = newFaceMapping[fid];
//        }
//        // =================== 重新映射结束 ===================
//    }
//
//    // 根据 cell 中包含的面信息，将单元划分为内部单元或边界单元
//    void ClassifySolidMatrixCells() 
//    {
//        for (const auto& cell : cells) 
//        {
//			bool isBoundary = false;// 边界单元标志
//            for (int faceId : cell.faceIDs) 
//            {
//                if (faces[faceId - 1].neighborCell == -1) { // 注意：如果 faces 数组索引从0开始，则需要减1
//                    isBoundary = true;
//                    break;
//                }
//            }
//            if (isBoundary)
//                boundaryCells.push_back(cell);
//            else
//                innerCells.push_back(cell);
//        }
//    }
//
//    // 输出网格基本信息
//    void printMeshInfo() {
//        cout << "总节点数: " << nodes.size() << endl;
//        cout << "总单元数（三角形）: " << cells.size() << endl;
//        cout << "总面数: " << faces.size() << endl;
//        cout << "内部单元数: " << innerCells.size() << endl;
//        cout << "边界单元数: " << boundaryCells.size() << endl;
//    }
//};
