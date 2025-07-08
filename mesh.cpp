#include "Mesh.h"
#include <iostream>
#include <algorithm>
#include <map>
#include <gmsh.h>
#include <set>
#include <cmath>
#include "Bound.h"

Mesh::Mesh() : gridCount_(0) {}


void Mesh::BuildMesh(double lengthX, double lengthY, double lengthZ, int nx, int ny, int nz, bool usePrism, bool useQuadBase)
{
    gmsh::initialize();
    gmsh::model::add("UnstructuredMesh-EDFM");

    bool is2D = (lengthZ <= 0.0);
    double lc = -1;

    if (is2D)
    {
        gmsh::option::setNumber("Mesh.RecombineAll", 0);
        gmsh::option::setNumber("Mesh.Algorithm", 6);
        lc = min(lengthX / nx, lengthY / ny);
        int p1 = gmsh::model::geo::addPoint(0, 0, 0, lc);
        int p2 = gmsh::model::geo::addPoint(lengthX, 0, 0, lc);
        int p3 = gmsh::model::geo::addPoint(lengthX, lengthY, 0, lc);
        int p4 = gmsh::model::geo::addPoint(0, lengthY, 0, lc);

        int l1 = gmsh::model::geo::addLine(p1, p2);
        int l2 = gmsh::model::geo::addLine(p2, p3);
        int l3 = gmsh::model::geo::addLine(p3, p4);
        int l4 = gmsh::model::geo::addLine(p4, p1);

        int loop = gmsh::model::geo::addCurveLoop({ l1, l2, l3, l4 });
        gmsh::model::geo::addPlaneSurface({ loop });
    }
    else if (usePrism)   ////存在问题，仍待解决
    {
        
        lc = min({ lengthX / nx, lengthY / ny, lengthZ / nz });
        // 关键设置：关闭自动细分，避免棱柱被拆成四面体
        gmsh::option::setNumber("Mesh.RecombineAll", 0);
        gmsh::option::setNumber("Mesh.SubdivisionAlgorithm", 0);
        gmsh::option::setNumber("Mesh.ElementOrder", 1); // 线性单元
        // 定义底面矩形四个点
        int p1 = gmsh::model::geo::addPoint(0, 0, 0, lc);
        int p2 = gmsh::model::geo::addPoint(lengthX, 0, 0, lc);
        int p3 = gmsh::model::geo::addPoint(lengthX, lengthY, 0, lc);
        int p4 = gmsh::model::geo::addPoint(0, lengthY, 0, lc);
        // 定义四条边
        int l1 = gmsh::model::geo::addLine(p1, p2);
        int l2 = gmsh::model::geo::addLine(p2, p3);
        int l3 = gmsh::model::geo::addLine(p3, p4);
        int l4 = gmsh::model::geo::addLine(p4, p1);
        int loop = gmsh::model::geo::addCurveLoop({ l1, l2, l3, l4 });
        int surface = gmsh::model::geo::addPlaneSurface({ loop });
        // 设置结构化网格分段数
        gmsh::model::geo::mesh::setTransfiniteCurve(l1, nx + 1);
        gmsh::model::geo::mesh::setTransfiniteCurve(l3, nx + 1);
        gmsh::model::geo::mesh::setTransfiniteCurve(l2, ny + 1);
        gmsh::model::geo::mesh::setTransfiniteCurve(l4, ny + 1);
        gmsh::model::geo::mesh::setTransfiniteSurface(surface);
        gmsh::model::geo::mesh::setTransfiniteSurface(surface, "Left", { p1, p2, p3, p4 });
        // 是否将底面四边形转成三角形（生成棱柱或六面体）
        if (useQuadBase)
            gmsh::model::geo::mesh::setRecombine(0, surface);
        // 执行挤压，最后一个参数控制是否生成六面体（true）或棱柱（false）
        gmsh::vectorpair surfaceVec = { {2, surface} };
        gmsh::vectorpair outDimTags;
        std::vector<int> numElements = { nz };
        gmsh::model::geo::extrude(surfaceVec, 0, 0, lengthZ, outDimTags, numElements, {}, /*recombine=*/useQuadBase);
    }
    else
    {
        gmsh::option::setNumber("Mesh.RecombineAll", 0);
        gmsh::option::setNumber("Mesh.Algorithm3D", 1);
        lc = std::min({ lengthX / nx, lengthY / ny, lengthZ / nz });

        int p1 = gmsh::model::geo::addPoint(0, 0, 0, lc);
        int p2 = gmsh::model::geo::addPoint(lengthX, 0, 0, lc);
        int p3 = gmsh::model::geo::addPoint(lengthX, lengthY, 0, lc);
        int p4 = gmsh::model::geo::addPoint(0, lengthY, 0, lc);
        int p5 = gmsh::model::geo::addPoint(0, 0, lengthZ, lc);
        int p6 = gmsh::model::geo::addPoint(lengthX, 0, lengthZ, lc);
        int p7 = gmsh::model::geo::addPoint(lengthX, lengthY, lengthZ, lc);
        int p8 = gmsh::model::geo::addPoint(0, lengthY, lengthZ, lc);

        int l1 = gmsh::model::geo::addLine(p1, p2);
        int l2 = gmsh::model::geo::addLine(p2, p3);
        int l3 = gmsh::model::geo::addLine(p3, p4);
        int l4 = gmsh::model::geo::addLine(p4, p1);

        int l5 = gmsh::model::geo::addLine(p5, p6);
        int l6 = gmsh::model::geo::addLine(p6, p7);
        int l7 = gmsh::model::geo::addLine(p7, p8);
        int l8 = gmsh::model::geo::addLine(p8, p5);

        int l9 = gmsh::model::geo::addLine(p1, p5);
        int l10 = gmsh::model::geo::addLine(p2, p6);
        int l11 = gmsh::model::geo::addLine(p3, p7);
        int l12 = gmsh::model::geo::addLine(p4, p8);

        int s1 = gmsh::model::geo::addCurveLoop({ l1, l2, l3, l4 });
        int s2 = gmsh::model::geo::addCurveLoop({ l5, l6, l7, l8 });
        int s3 = gmsh::model::geo::addCurveLoop({ l1, l10, -l5, -l9 });
        int s4 = gmsh::model::geo::addCurveLoop({ l2, l11, -l6, -l10 });
        int s5 = gmsh::model::geo::addCurveLoop({ l3, l12, -l7, -l11 });
        int s6 = gmsh::model::geo::addCurveLoop({ l4, l9, -l8, -l12 });

        int sf1 = gmsh::model::geo::addPlaneSurface({ s1 });
        int sf2 = gmsh::model::geo::addPlaneSurface({ s2 });
        int sf3 = gmsh::model::geo::addPlaneSurface({ s3 });
        int sf4 = gmsh::model::geo::addPlaneSurface({ s4 });
        int sf5 = gmsh::model::geo::addPlaneSurface({ s5 });
        int sf6 = gmsh::model::geo::addPlaneSurface({ s6 });

        int sl = gmsh::model::geo::addSurfaceLoop({ sf1, sf2, sf3, sf4, sf5, sf6 });
        gmsh::model::geo::addVolume({ sl });
    }

    gmsh::model::geo::synchronize();
    gmsh::model::mesh::generate(is2D ? 2 : 3);
	


    // === Step 1: 获取所有节点，并建立 ID → 坐标 的映射 ===
    vector<size_t> nodeTags;
    vector<double> nodeCoords, parametricCoords;
    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, parametricCoords);

    for (size_t i = 0; i < nodeTags.size(); i++)
    {
		int id = static_cast<int>(nodeTags[i]); // 获取节点编号
		cout << "Node ID = " << id << endl;
        Vector coord(nodeCoords[3 * i], nodeCoords[3 * i + 1], nodeCoords[3 * i + 2]);
        nodes_.emplace_back(id, coord);        // 括号法定义，添加到节点列表
        nodesMap_[id] = Node(id, coord);        // 建立映射关系
    }

    // ==== Step 2: 读取单元信息并构建 Cell（支持三角形和四面体） ====
	vector<int> elementTypes;  //储存单元类型
	vector<vector<size_t>> elementTags, nodeTagsPerElement;  //二维数组，按照单元类型储存单元标签和组成单元的节点信息
	gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTagsPerElement, is2D ? 2 : 3);  //is2D-True-2D网格，False-3D网格

    for (size_t i = 0; i < elementTypes.size(); i++)
    {
        int type = elementTypes[i]; //获取网格类型
        cout << "Element type = " << elementTypes[i] << endl;

        // 三角形单元（2D） 1=线单元, 2=三角形, 3=四边形  其他网格类型待补充
        if (is2D && type == 2)
        {
            size_t numTriangles = elementTags[i].size();
            cout << "Element numbers = " << numTriangles << endl;
            for (size_t j = 0; j < numTriangles; j++)               //遍历每一个网格单元
            {
				int cellId = static_cast<int>(elementTags[i][j]);   // 读取单元编号
				cout << "Cell ID = " << cellId << endl;
				vector<int> cellNodeIDs;                            //建立储存单元节点编号的容器
                for (int k = 0; k < 3; k++)                         // 三角形单元有三个节点
                cellNodeIDs.push_back(static_cast<int>(nodeTagsPerElement[i][3 * j + k]));               
                Cell cell(cellId, cellNodeIDs);
                cell.computeCenterAndVolume(nodesMap_);
                cells_.push_back(cell);
            }
        }

        // 四面体单元（3D）4=四面体 5=六面体 6=棱柱 7=金字塔
        else if (!is2D && type == 4)
        {
            size_t numTets = elementTags[i].size();
            for (size_t j = 0; j < numTets; j++)
            {
				int cellId = static_cast<int>(elementTags[i][j]); // 读取单元编号
				vector<int> cellNodeIDs;                         //建立储存单元节点编号的容器
				for (int k = 0; k < 4; k++)                      // 四面体单元有四个节点
                    cellNodeIDs.push_back(static_cast<int>(nodeTagsPerElement[i][4 * j + k]));

				Cell cell(cellId, cellNodeIDs);       // 构造 Cell 对象
				cell.computeCenterAndVolume(nodesMap_); // 计算 Cell 的中心和体积
                cells_.push_back(cell);
            }
        }
        // 棱柱单元（type == 6）
        else if (!is2D && type == 6)
        {
            size_t numPrisms = elementTags[i].size();
            for (size_t j = 0; j < numPrisms; j++)
            {
                int cellId = static_cast<int>(elementTags[i][j]);
                vector<int> cellNodeIDs;
                for (int k = 0; k < 6; k++)  // 棱柱有 6 个节点
                    cellNodeIDs.push_back(static_cast<int>(nodeTagsPerElement[i][6 * j + k]));

                Cell cell(cellId, cellNodeIDs);
                cell.computeCenterAndVolume(nodesMap_);
                cells_.push_back(cell);
            }
        }

        // 六面体单元（type == 5）
        else if (!is2D && type == 5)
        {
            size_t numPrisms = elementTags[i].size();
            for (size_t j = 0; j < numPrisms; j++)
            {
                int cellId = static_cast<int>(elementTags[i][j]);
                vector<int> cellNodeIDs;
                for (int k = 0; k < 8; k++)  // 棱柱有 6 个节点
                    cellNodeIDs.push_back(static_cast<int>(nodeTagsPerElement[i][8 * j + k]));

                Cell cell(cellId, cellNodeIDs);
                cell.computeCenterAndVolume(nodesMap_);
                cells_.push_back(cell);
            }
        }


    }
	//构建 id → 下标 映射，方便后续快速访问
    cellId2index_.clear();
    for (size_t i = 0; i < cells_.size(); ++i) {

        cout << "Cell ID: " << cells_[i].id << endl;
			
        cellId2index_[cells_[i].id] = static_cast<int>(i);
    }
    
    gridCount_ = cells_.size();  // 统计总单元数

	// === Step 3: 构造面（Face）并建立 Cell 和 Face 的对应关系 ===  *注意 对于3D情况而言，是构造面片与单元的对应关系

	//map<set<int>, vector<int>> faceMap;   //set的作用是按从小到大自动排序，去除重复情况，可以当作唯一标识某个面的 key  O(logn)  后续可用HashMap优化
 //   for (const auto& cell : cells_)
 //   {
	//	const auto& cn = cell.nodeIDs; // 获取单元的节点编号
 //      
 //       if (cn.size() == 3)
 //       {
 //           // 2D 情况（三角形边）
 //           vector<set<int>> edges =
 //           {
 //               {cn[0], cn[1]},
 //               {cn[1], cn[2]},
 //               {cn[2], cn[0]}
 //           };
 //           for (const auto& edge : edges)
 //               faceMap[edge].push_back(cell.id);
 //       }
 //      
 //       else if (cn.size() == 4)
 //       {
 //           // 3D情况四面体的四个面（每面是三个点）
 //           vector<set<int>> faces = 
 //           {
 //               {cn[0], cn[1], cn[2]},
 //               {cn[0], cn[1], cn[3]},
 //               {cn[1], cn[2], cn[3]},
 //               {cn[0], cn[2], cn[3]}
 //           };
 //           for (const auto& face : faces)
 //               faceMap[face].push_back(cell.id);
 //       }

 //       else if (cn.size() == 6)  // 棱柱
 //       {
 //           vector<set<int>> faces =
 //           {
 //               {cn[0], cn[1], cn[2]}, // 三角底面
 //               {cn[3], cn[4], cn[5]}, // 三角顶面
 //               {cn[0], cn[1], cn[4], cn[3]}, // 四边形面
 //               {cn[1], cn[2], cn[5], cn[4]}, // 四边形面
 //               {cn[2], cn[0], cn[3], cn[5]}  // 四边形面
 //           };
 //           for (const auto& face : faces)
 //           faceMap[face].push_back(cell.id);
 //              
 //           
 //       }

 //       else if (cn.size() == 8)  // 六面体 Hexahedron
 //       {
 //           std::vector<std::set<int>> faces =
 //           {
 //               {cn[0], cn[1], cn[2], cn[3]}, // 底面
 //               {cn[4], cn[5], cn[6], cn[7]}, // 顶面
 //               {cn[0], cn[1], cn[5], cn[4]}, // 侧面1
 //               {cn[1], cn[2], cn[6], cn[5]}, // 侧面2
 //               {cn[2], cn[3], cn[7], cn[6]}, // 侧面3
 //               {cn[3], cn[0], cn[4], cn[7]}  // 侧面4
 //           };
 //           for (const auto& face : faces)
 //               faceMap[face].push_back(cell.id);
 //       }      
 //   }

	//// 面的编号接着Cell的编号继续
 //   int faceId = 1;

	//// 遍历面映射，构造 Face 对象并建立 Cell 和 Face 的对应关系
 //   for (const auto& entry : faceMap) 
 //   {
	//	vector<int> nodeIDs(entry.first.begin(), entry.first.end()); //遍历map并存储面节点编号
 //       vector<Vector> nodeCoords;
 //       for (int nid : nodeIDs)
 //           nodeCoords.push_back(nodesMap_[nid].coord);

 //       Face face(faceId, nodeIDs, nodeCoords);
	//	face.ownerCell = entry.second[0];   //第一个储存的单元编号为拥有该面的单元编号
	//	face.neighborCell = (entry.second.size() == 2) ? entry.second[1] : -1; //如果是边界面那么他的neighborCell为-1
 //       faces_.push_back(face);

 //       for (int cid : entry.second)
 //       {
 //           for (auto& cell : cells_) 
 //           {
 //               if (cell.id == cid) 
 //               {
 //                   cell.faceIDs.push_back(faceId);
 //                   break;
 //               }
 //           }
 //       }
 //       faceId++;
 //   }


	// === Step 4: 构造 Face 对象并建立 Cell 和 Face 的对应关系 ===
    std::unordered_map<std::vector<int>, std::vector<int>, VectorHash> faceMap;

    for (const auto& cell : cells_) {
        for (const auto& faceNodeIDs : cell.getLocalFaces()) {
            std::vector<int> sortedNodeIDs = faceNodeIDs;
            std::sort(sortedNodeIDs.begin(), sortedNodeIDs.end());
            faceMap[sortedNodeIDs].push_back(cell.id);
        }
    }

    faces_.clear();
    int faceId = 1;
    for (const auto& entry : faceMap)
    {
        const std::vector<int>& nodeIDs = entry.first;
        std::vector<Vector> nodeCoords;
        for (int nid : nodeIDs)
            nodeCoords.push_back(nodesMap_[nid].coord);

        Face face(faceId, nodeIDs, nodeCoords);
        face.ownerCell = entry.second[0];
        face.neighborCell = (entry.second.size() == 2) ? entry.second[1] : -1;
        faces_.push_back(face);

        for (int cid : entry.second)
        {
            int idx = cellId2index_.at(cid);
            cells_[idx].faceIDs.push_back(faceId);
        }
        faceId++;
    }
    std::cout << "BuildMesh: faces_ count = " << faces_.size() << std::endl;
    std::cout << "=== Face owner/neighbor check ===" << std::endl;
    for (const auto& face : faces_) {
        std::cout << "Face " << face.id << ": [";
        for (int nid : face.nodeIDs)
            std::cout << nid << " ";
        std::cout << "], owner = " << face.ownerCell
            << ", neighbor = " << face.neighborCell << std::endl;
    }
    gmsh::write("UnstructuredMesh-EDFM.msh");
    gmsh::finalize();
}


void Mesh::ClassifySolidMatrixCells()
{
    // 构造一个 faceId 到 Face 的映射
   
    for (auto& cell : cells_)
    {
        bool isBoundary = false;
        for (int faceId : cell.faceIDs)
        {
            if (faces_[faceId - 1].neighborCell == -1)
            {
                isBoundary = true;
                break;
            }
        }
        cell.location = isBoundary ? Cell::LocationType::Boundary : Cell::LocationType::Inner;
    }
}

std::vector<std::reference_wrapper<const Cell>> Mesh::getInnerCells() const
{
    std::vector<std::reference_wrapper<const Cell>> innerCells;
    for (const auto& cell : cells_)
    {
        if (cell.location == Cell::LocationType::Inner)
            innerCells.push_back(std::cref(cell));
    }
    return innerCells;
}

std::vector<std::reference_wrapper<const Cell>> Mesh::getBoundaryCells() const
{
	std::vector<std::reference_wrapper<const Cell>> boundaryCells;
	for (const auto& cell : cells_)
	{
		if (cell.location == Cell::LocationType::Boundary)
			boundaryCells.push_back(std::cref(cell));
	}
	return boundaryCells;
}

void Mesh::CreateSolidMatrixGhostCells() 
{
    int origN = (int)cells_.size();
    ghostStartIndex_ = origN;  // 0-based

    std::vector<Cell> ghosts;
    ghosts.reserve(faces_.size());

    int ghostCount = 0;
    // 1) 遍历所有面，给每个边界面造一个 ghost cell，并更新 face.neighborCell
    for (auto& face : faces_) 
    {
        if (!face.isBoundary()) continue;

        // 拿出它的 ownerCell 在 cells 中的下标
        int ownerID = face.ownerCell;
		int idx = cellId2index_.at(ownerID); //找到 cellId 对应的下标
        Cell   realC = cells_[idx];

        // 构造 ghost
        Cell ghostC = realC;
        ++ghostCount;
        // 给它一个新 ID（不能跟已有 ID 冲突）
        ghostC.id = -ghostCount;
        // 更新那条 Face 的 neighborCell
        face.neighborCell = ghostC.id;

        // === （你的镜像逻辑，把 ghostC.center 写到面外） ===
        Vector Cp = realC.center;
        Vector mp = face.midpoint;
        Vector n = face.normal;      // 单位法向
        Vector v = Cp - mp;
        double d = v * n;
        ghostC.center = Cp - 2.0 * d * n;

        // ghost 不需要原来的 nodeIDs/faceIDs
        ghostC.nodeIDs.clear();
        ghostC.faceIDs.clear();
        ghosts.push_back(ghostC);
    }

    // 2) 把 ghosts append 到 cells，并更新 cellId2index
    for (auto& g : ghosts) 
    {
        cellId2index_[g.id] = (int)cells_.size();
        cells_.push_back(g);
    }
}

void Mesh::printMeshInfo()
{
    cout << "总节点数: " << nodes_.size() << endl;
    cout << "总单元数: " << cells_.size() << endl;
    cout << "总面数: " << faces_.size() << endl;
    int nInner = 0, nBoundary = 0;
    for (const auto& cell : cells_) {
        if (cell.location == Cell::LocationType::Inner) ++nInner;
        else if (cell.location == Cell::LocationType::Boundary) ++nBoundary;
    }
    cout << "内部单元数: " << nInner << endl;
    cout << "边界单元数: " << nBoundary << endl;

    cout << "\n―― 节点信息 ――" << endl;
    for (const auto& node : nodes_)
    {
        cout << "Node " << node.id << " : ("
            << node.coord.m_x << ", " << node.coord.m_y << ", " << node.coord.m_z << ")"
            << endl;
    }
    cout << "\n―― 单元信息 ――" << endl;      //还需要添加Owner和Neighbor
    for (const auto& cell : cells_)
    {
        if (cell.id < 0) continue;   // 跳过ghost cell
        cout << "Cell " << cell.id << " 中心: ("
            << cell.center.m_x << ", " << cell.center.m_y << ", " << cell.center.m_z << ")"
            << " 面积: " << cell.volume << endl;
        cout << "  包含面: ";
        for (int fid : cell.faceIDs)
            cout << fid << " ";
        cout << endl;
    }
    cout << "\n―― 面信息 ――" << endl;
    for (const auto& face :faces_)
    {
        cout << "Face " << face.id << " 节点: ";
        for (int nid : face.nodeIDs)
            cout << nid << " ";
        cout << " | owner: " << face.ownerCell
            << " neighbor: " << face.neighborCell << endl;
        cout << "  中点: (" << face.midpoint.m_x << ", " << face.midpoint.m_y << ", " << face.midpoint.m_z
            << "), 边长: " << face.length << endl;
    }

}

double Mesh::getCellArea(int cellID) const
{
    // 根据 cellID 查找基岩单元，获取它的面积
    for (const auto& cell : cells_)
    {
        if (cell.id == cellID) 
        {
            return cell.volume;  
        }
    }
    return 0.0; // 如果找不到对应的单元，返回0
}

void Mesh::ComputeSolidMatrixMeshFaceGeometricInfor(NormalVectorCorrectionMethod method)
{
    // 遍历所有面
    for (auto& face : faces_) 
    {
        // 只处理内部面
        //if (face.isBoundary()) continue;

        // 找到 owner 和 neighbor 在 cells 向量中的下标
        auto it_o = cellId2index_.find(face.ownerCell);
        auto it_n = cellId2index_.find(face.neighborCell);
        if (it_o == cellId2index_.end() || it_n == cellId2index_.end()) 
        {
            // 安全起见，跳过或报错

            continue;
        }

        // 提取两个单元的几何中心
        const Vector& Cp = cells_[it_o->second].center;
        const Vector& Cn = cells_[it_n->second].center;

        // 调用分解函数，能选三种方法之一
        face.computeFaceVectors(Cp, Cn, method);
    }
}

// 导出 ghostStartIndex 到一个单独的 txt
void Mesh::exportGhostInfo(const std::string& prefix) const 
{
    std::ofstream f(prefix + "_info.txt");
    // 仅写一个数字：第一个 ghost cell 在 cells[] 中的下标
    f << ghostStartIndex_ << "\n";
    f.close();
}

void Mesh::exportToTxt(const std::string& prefix) const 
{
    // Nodes
    ofstream fn(prefix + "_nodes.txt");
    fn << "id x y z\n";
    for (auto& n : nodes_)
        fn << n.id << " "
        << n.coord.m_x << " "
        << n.coord.m_y << " "
        << n.coord.m_z << "\n";
    fn.close();

    // Faces
    ofstream ff(prefix + "_faces.txt");
    ff << "id n1 n2 mx my mz\n";
    for (auto& f : faces_) {
        ff << f.id << " "
            << f.nodeIDs[0] << " "
            << f.nodeIDs[1] << " "
            << f.midpoint.m_x << " "
            << f.midpoint.m_y << " "
            << f.midpoint.m_z << "\n";
    }
    ff.close();

    // Cells
    ofstream fc(prefix + "_cells.txt");
    fc << "id cx cy cz\n";
    for (auto& c : cells_) {
        fc << c.id << " "
            << c.center.m_x << " "
            << c.center.m_y << " "
            << c.center.m_z << "\n";
    }
    fc.close();

    exportGhostInfo(prefix);
}
