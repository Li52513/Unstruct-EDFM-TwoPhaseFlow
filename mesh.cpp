#include "Mesh.h"
#include <iostream>
#include <algorithm>
#include <map>
#include <gmsh.h>
#include <set>
#include <cmath>
#include "Bound.h"
#include <unordered_set>

Mesh::Mesh() : gridCount_(0) {}



void Mesh::BuildMesh(double lengthX, double lengthY, double lengthZ, int nx, int ny, int nz, bool usePrism, bool useQuadBase)
{
    //// 初始化GMSH环境
    //gmsh::initialize();
    //gmsh::model::add("UnstructuredMesh-EDFM");
    //// 判断是否为2D网格（Z方向长度≤0）
    //bool is2D = (lengthZ <= 0.0);
    //double lc = -1; ;  // 特征长度，初始化为无效值
    //// ----网格剖分：构造一个矩形面（2D面 3D为底面）
    //if (is2D)
    //{
    //    // 2D网格设置：禁用网格重组，使用网格算法6
    //    gmsh::option::setNumber("Mesh.RecombineAll", 0);  //三角形网格 0 四边形网格 1
    //    gmsh::option::setNumber("Mesh.Algorithm", 8);   //三角形网格 6 1-8
    //    // 计算特征长度（取X和Y方向的最小网格尺寸）
    //    lc = min(lengthX / nx, lengthY / ny);
    //    // 创建矩形四个顶点
    //    int p1 = gmsh::model::geo::addPoint(0, 0, 0, lc);
    //    int p2 = gmsh::model::geo::addPoint(lengthX, 0, 0, lc);
    //    int p3 = gmsh::model::geo::addPoint(lengthX, lengthY, 0, lc);
    //    int p4 = gmsh::model::geo::addPoint(0, lengthY, 0, lc);
    //    // 创建四条边
    //    int l1 = gmsh::model::geo::addLine(p1, p2);
    //    int l2 = gmsh::model::geo::addLine(p2, p3);
    //    int l3 = gmsh::model::geo::addLine(p3, p4);
    //    int l4 = gmsh::model::geo::addLine(p4, p1);
    //    // 创建曲线环和平面表面
    //    int loop = gmsh::model::geo::addCurveLoop({ l1, l2, l3, l4 });
    //    gmsh::model::geo::addPlaneSurface({ loop });
    //}
    //else if (usePrism)   ///< 棱柱单元生成（实验性功能，存在问题）
    //{
    //    // 计算特征长度（取X、Y和Z方向的最小网格尺寸）
    //    lc = min({ lengthX / nx, lengthY / ny, lengthZ / nz });
    //    // 关键设置：关闭自动细分，避免棱柱被拆成四面体
    //    gmsh::option::setNumber("Mesh.RecombineAll", 0);
    //    gmsh::option::setNumber("Mesh.SubdivisionAlgorithm", 0);
    //    gmsh::option::setNumber("Mesh.ElementOrder", 1); // 线性单元
    //    // 定义底面矩形四个点
    //    int p1 = gmsh::model::geo::addPoint(0, 0, 0, lc);
    //    int p2 = gmsh::model::geo::addPoint(lengthX, 0, 0, lc);
    //    int p3 = gmsh::model::geo::addPoint(lengthX, lengthY, 0, lc);
    //    int p4 = gmsh::model::geo::addPoint(0, lengthY, 0, lc);
    //    // 定义四条边
    //    int l1 = gmsh::model::geo::addLine(p1, p2);
    //    int l2 = gmsh::model::geo::addLine(p2, p3);
    //    int l3 = gmsh::model::geo::addLine(p3, p4);
    //    int l4 = gmsh::model::geo::addLine(p4, p1);
    //    int loop = gmsh::model::geo::addCurveLoop({ l1, l2, l3, l4 });
    //    int surface = gmsh::model::geo::addPlaneSurface({ loop });
    //   // 设置结构化网格分段数
    //    gmsh::model::geo::mesh::setTransfiniteCurve(l1, nx + 1);
    //    gmsh::model::geo::mesh::setTransfiniteCurve(l3, nx + 1);
    //    gmsh::model::geo::mesh::setTransfiniteCurve(l2, ny + 1);
    //    gmsh::model::geo::mesh::setTransfiniteCurve(l4, ny + 1);
    //    gmsh::model::geo::mesh::setTransfiniteSurface(surface);
    //    gmsh::model::geo::mesh::setTransfiniteSurface(surface, "Left", { p1, p2, p3, p4 });
    //    // 是否将底面四边形转成三角形（生成棱柱或六面体）
    //    if (useQuadBase)
    //        gmsh::model::geo::mesh::setRecombine(1, surface);
    //    // 执行挤压，最后一个参数控制是否生成六面体（true）或棱柱（false）
    //    gmsh::vectorpair surfaceVec = { {2, surface} };
    //    gmsh::vectorpair outDimTags;
    //    std::vector<int> numElements = { nz };
    //    gmsh::model::geo::extrude(surfaceVec, 0, 0, lengthZ, outDimTags, numElements, {}, /*recombine=*/useQuadBase);
    //}
    //else
    //{
    //    // 3D四面体网格设置
    //    gmsh::option::setNumber("Mesh.RecombineAll", 0);
    //    gmsh::option::setNumber("Mesh.Algorithm3D", 1);
    //    // 计算特征长度
    //    lc = std::min({ lengthX / nx, lengthY / ny, lengthZ / nz });
    //    // 创建立方体八个顶点
    //    int p1 = gmsh::model::geo::addPoint(0, 0, 0, lc);
    //    int p2 = gmsh::model::geo::addPoint(lengthX, 0, 0, lc);
    //    int p3 = gmsh::model::geo::addPoint(lengthX, lengthY, 0, lc);
    //    int p4 = gmsh::model::geo::addPoint(0, lengthY, 0, lc);
    //    int p5 = gmsh::model::geo::addPoint(0, 0, lengthZ, lc);
    //    int p6 = gmsh::model::geo::addPoint(lengthX, 0, lengthZ, lc);
    //    int p7 = gmsh::model::geo::addPoint(lengthX, lengthY, lengthZ, lc);
    //    int p8 = gmsh::model::geo::addPoint(0, lengthY, lengthZ, lc);
    //    // 创建十二条边
    //    int l1 = gmsh::model::geo::addLine(p1, p2);
    //    int l2 = gmsh::model::geo::addLine(p2, p3);
    //    int l3 = gmsh::model::geo::addLine(p3, p4);
    //    int l4 = gmsh::model::geo::addLine(p4, p1);
    //    int l5 = gmsh::model::geo::addLine(p5, p6);
    //    int l6 = gmsh::model::geo::addLine(p6, p7);
    //    int l7 = gmsh::model::geo::addLine(p7, p8);
    //    int l8 = gmsh::model::geo::addLine(p8, p5);
    //    int l9 = gmsh::model::geo::addLine(p1, p5);
    //    int l10 = gmsh::model::geo::addLine(p2, p6);
    //    int l11 = gmsh::model::geo::addLine(p3, p7);
    //    int l12 = gmsh::model::geo::addLine(p4, p8);
    //    // 创建六个表面
    //    int s1 = gmsh::model::geo::addCurveLoop({ l1, l2, l3, l4 });
    //    int s2 = gmsh::model::geo::addCurveLoop({ l5, l6, l7, l8 });
    //    int s3 = gmsh::model::geo::addCurveLoop({ l1, l10, -l5, -l9 });
    //    int s4 = gmsh::model::geo::addCurveLoop({ l2, l11, -l6, -l10 });
    //    int s5 = gmsh::model::geo::addCurveLoop({ l3, l12, -l7, -l11 });
    //    int s6 = gmsh::model::geo::addCurveLoop({ l4, l9, -l8, -l12 });
    //    // 创建平面表面
    //    int sf1 = gmsh::model::geo::addPlaneSurface({ s1 });
    //    int sf2 = gmsh::model::geo::addPlaneSurface({ s2 });
    //    int sf3 = gmsh::model::geo::addPlaneSurface({ s3 });
    //    int sf4 = gmsh::model::geo::addPlaneSurface({ s4 });
    //    int sf5 = gmsh::model::geo::addPlaneSurface({ s5 });
    //    int sf6 = gmsh::model::geo::addPlaneSurface({ s6 });
    //    // 创建表面环和体积
    //    int sl = gmsh::model::geo::addSurfaceLoop({ sf1, sf2, sf3, sf4, sf5, sf6 });
    //    gmsh::model::geo::addVolume({ sl });
    //}
    //// 同步几何模型并生成网格
    //gmsh::model::geo::synchronize();
    //gmsh::model::mesh::generate(is2D ? 2 : 3);
    
    
	// 0) 初始化 Gmsh
    gmsh::initialize();
    gmsh::model::add("UnstructuredMesh-EDFM");

    // 视作 2D 的条件：长度Z<=0 或者 nz==0
    const bool is2D = (lengthZ <= 0.0 || nz <= 0);
    double lc = -1.0;

    // ---- 通用：构造一个矩形底面几何（2D面 or 3D的底面） ----
    {
        // 尺寸：按目标网格数估算平均特征长度
        lc = std::min(lengthX / std::max(nx, 1), lengthY / std::max(ny, 1));
        // 顶点
        int p1 = gmsh::model::geo::addPoint(0, 0, 0, lc);
        int p2 = gmsh::model::geo::addPoint(lengthX, 0, 0, lc);
        int p3 = gmsh::model::geo::addPoint(lengthX, lengthY, 0, lc);
        int p4 = gmsh::model::geo::addPoint(0, lengthY, 0, lc);

        // 边
        int l1 = gmsh::model::geo::addLine(p1, p2);
        int l2 = gmsh::model::geo::addLine(p2, p3);
        int l3 = gmsh::model::geo::addLine(p3, p4);
        int l4 = gmsh::model::geo::addLine(p4, p1);

        // 面
        int loop = gmsh::model::geo::addCurveLoop({ l1, l2, l3, l4 });
        int surface = gmsh::model::geo::addPlaneSurface({ loop });

        // ---- 二维网格类型：三角形 or 四边形（非结构化） ----
        if (useQuadBase)
        {
            //非结构化四边形
            gmsh::option::setNumber("Mesh.RecombineAll", 1);  // 三角→四边形重组
            gmsh::option::setNumber("Mesh.Algorithm", 4);     // 8Quad 前沿法（更稳）
            gmsh::option::setNumber("Mesh.RecombinationAlgorithm", 2); // blossom
        }
        else
        {
            // 非结构化三角
            gmsh::option::setNumber("Mesh.RecombineAll", 0);
            gmsh::option::setNumber("Mesh.Algorithm", 6);     // Frontal-Delaunay
        }
        // 优化
        gmsh::option::setNumber("Mesh.Optimize", 1);
        gmsh::option::setNumber("Mesh.Smoothing", 2);

        //改动：添加
        gmsh::model::geo::synchronize();
        gmsh::model::mesh::generate(2);

        //------如果需要生成三维网格，则沿z方向扫掠---
        if (!is2D && usePrism)
        {
            //改动：添加
            gmsh::option::setNumber("Mesh.SubdivisionAlgorithm", 0);
            gmsh::option::setNumber("Mesh.ElementOrder", 1);

            gmsh::vectorpair base = { {2, surface} };
            gmsh::vectorpair outDimTags;

            // 非均匀层厚：如果 setLayerHeights() 提供了配置(待添加），就用它；否则均分 nz
            std::vector<int>    numElems;   // 通常给 {nz}
            std::vector<double> heights;    // 或者给每层厚度
            numElems = { std::max(nz, 1) };  //均分z方向网格

            // 选择是否重组：三角底 → 棱柱（recombine=false）；四边形底 → 六面体网格（recombine=true），由于gmsh在拉伸生成非结构化棱柱网格时候底面只能为三角形网格
            bool recombineForHex = useQuadBase;
            if (!useQuadBase) {
                std::cout << "[INFO] 3D 扫掠：三角底 → 生成棱柱(type=6)。" << std::endl;
            }
            else {
                std::cout << "[INFO] 3D 扫掠：四边形底 → 生成六面体(type=5)。" << std::endl;
            }
           
            //扫掠
            gmsh::model::geo::extrude
            (
                base, 0, 0, lengthZ,
                outDimTags,
                numElems,            // 均分层数
                heights,             // 或者非均匀层厚
                /*recombine=*/true   //将三角形进行合并 才能生成棱柱网格，如果是false 当底面为非结构化四边形的时候也可以生成六面体网格，但是当底面为三角形时生成的史四面体网格。 三角底: true→棱柱(type=6), false→四面体(type=4) 四边底: true→六面体(type=5), false→四面体(type=4)
            );

            gmsh::model::geo::synchronize();
            gmsh::model::mesh::generate(3);

        }

    }
    

    // === Step 1: 获取所有节点，并建立 ID → 坐标 的映射 ===

    vector<size_t> nodeTags;            ///< GMSH节点ID集合
    vector<double> nodeCoords;          ///< 节点坐标扁平化数组 [x1,y1,z1,x2,y2,z2,...]
    vector<double> parametricCoords;    ///< 参数坐标
    // 调用GMSH API获取网格节点数据
    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, parametricCoords);
    nodes_.clear(); nodes_.reserve(nodeTags.size());//清除nodes_并预留空间
    nodesMap_.clear(); nodesMap_.reserve(nodeTags.size()); //清除nodesMap_并预留空间

    // 处理所有节点数据
    for (size_t i = 0; i < nodeTags.size(); i++)
    {
        // 转换节点ID类型: size_t → int (GMSH → 本地系统)
		int id = static_cast<int>(nodeTags[i]); 

        // 提取节点三维坐标 (x, y, z)
		Vector coord(
            nodeCoords[3 * i],          // X坐标
            nodeCoords[3 * i + 1],      // Y坐标
			nodeCoords[3 * i + 2]		// Z坐标
        );

        // 将节点添加到线性容器（保持原始顺序）
        nodes_.emplace_back(id, coord);  

        // 将节点添加到映射表（便于通过ID快速查找）
        nodesMap_[id] = Node(id, coord);   
    }

    // ==== Step 2: 读取单元信息，识别组成单元的网格节点，计算终点坐标和体积大小，从而构建 Cell（支持三角形和四面体）数据类型====

	vector<int> elementTypes;                               //储存单元类型
	vector<vector<size_t>> elementTags, nodeTagsPerElement;  //二维数组，按照单元类型储存单元标签和组成单元的节点信息

    // 获取网格信息（is2D-True-2D网格，False-3D网格)
	gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTagsPerElement, is2D ? 2 : 3);  

    cells_.clear();
    // 处理所有单元类型
    for (size_t i = 0; i < elementTypes.size(); i++)
    {
        int type = elementTypes[i]; //获取网格类型
        std::cout << "Element type = " << elementTypes[i] << endl;

        // 三角形单元（2D） 1=线单元, 2=三角形, 3=四边形  其他网格类型待补充
        std::cout << "------------Cell全局编号和坐标信息------------" << endl;

        if (is2D && type == 2)
        { 
            size_t numTriangles = elementTags[i].size();    //有多少个三角形网格单元
            
            std::cout << "Element numbers = " << numTriangles << endl;
            
            //遍历每一个网格单元
            for (size_t j = 0; j < numTriangles; j++)               
            {
				int cellId = static_cast<int>(elementTags[i][j]);   // 读取网格单元局部编号
				//cout << "Cell ID = " << cellId << endl;
				vector<int> cellNodeIDs;                            //建立储存单元节点编号的容器  **建立cell与Node之间的连接关系
               
                // 提取三角形单元的三个节点
                for (int k = 0; k < 3; k++)
				{
                    // 三角形单元有三个节点，其中 k=0,1,2 分别对应三个节点 //局部编号
                
                    cellNodeIDs.push_back(static_cast<int>(nodeTagsPerElement[i][3 * j + k]));  
					cout << "Node ID = " << cellNodeIDs[k] << endl<<endl;
                    
				}
                
                // 创建Cell对象并计算几何属性
                Cell cell(cellId, cellNodeIDs);
                //cout <<"------编号为"<< cellId<<"的Cell中心点和面积信息--------"  << endl;
                cell.computeCenterAndVolume(nodesMap_);
                cells_.push_back(cell);
               // cout << "------------------------------------------------" << endl;
            }
        }
        else if (is2D && type == 3)
        {
            size_t n = elementTags[i].size();
            std::cout << "Element numbers = " << n << endl;
            for (size_t j = 0; j < n; ++j)
            {
                int cellId = static_cast<int>(elementTags[i][j]);
                vector<int> cellNodeIDs;
                for (int k = 0; k < 4; k++)
                {
                    // 四边形单元有四个节点，其中 k=0,1,2,3 分别对应四个节点 //局部编号

                    cellNodeIDs.push_back(static_cast<int>(nodeTagsPerElement[i][4 * j + k]));
                    cout << "Node ID = " << cellNodeIDs[k] << endl << endl;

                }
                // 创建Cell对象并计算几何属性
                Cell cell(cellId, cellNodeIDs);
                //cout <<"------编号为"<< cellId<<"的Cell中心点和面积信息--------"  << endl;
                cell.computeCenterAndVolume(nodesMap_);
                cells_.push_back(cell);
                // cout << "------------------------------------------------" << endl;
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

        else if (!is2D && type == 7)   // pyramid(5)（兼容读取）
        {         
            size_t n = elementTags[i].size();
            for (size_t j = 0; j < n; ++j) 
            {
                int cellId = static_cast<int>(elementTags[i][j]);
                std::vector<int> ids(5);
                for (int k = 0; k < 5; ++k) ids[k] = static_cast<int>(nodeTagsPerElement[i][5 * j + k]);
                Cell c(cellId, ids); c.computeCenterAndVolume(nodesMap_); cells_.push_back(c);
            }
        }


    }

    // 构建单元ID到索引的映射，方便后续快速访问
    cellId2index_.clear();
    for (size_t i = 0; i < cells_.size(); ++i)
    {
       // cout << "Cell ID: " << cells_[i].id << endl;	
		cellId2index_[cells_[i].id] = static_cast<int>(i); //局部编号到全局编号的映射
    }
    
    // 统计总单元数
    gridCount_ = cells_.size();
    std::cout << "网格数量为" << gridCount_ << endl;

    struct FaceAcc {
        std::vector<int> ordered;      // 有序顶点（来自某个单元的 getLocalFaces 原序）
        std::vector<int> cells;        // 邻接单元 id（1 或 2）
    };
    std::unordered_map<std::vector<int>, FaceAcc, VectorHash> FacetoCellMap;

    for (const auto& cell : cells_) {
        for (const auto& faceNodeIDs_ordered : cell.getLocalFaces(nodesMap_)) {
            std::vector<int> key = faceNodeIDs_ordered;
            std::sort(key.begin(), key.end());            // 只用于去重
            auto& acc = FacetoCellMap[key];
            if (acc.ordered.empty()) acc.ordered = faceNodeIDs_ordered; // 保存环序一次
            acc.cells.push_back(cell.id);
        }
    }

    // 2) 实例化：用“环序”构造几何；必要时校正朝向（可选）
    faces_.clear();
    int faceId = 1;
    for (auto& kv : FacetoCellMap) {
        const auto& orderedIDs = kv.second.ordered;   // 正确环序
        std::vector<Vector> nodeCoords; nodeCoords.reserve(orderedIDs.size());
        for (int nid : orderedIDs) nodeCoords.push_back(nodesMap_.at(nid).coord);

        Face face(faceId, orderedIDs, nodeCoords);    // 几何计算用环序
        face.ownerCell = kv.second.cells[0];
        face.neighborCell = (kv.second.cells.size() == 2) ? kv.second.cells[1] : -1;

        // （可选）法向一致化：让法向从 owner 指向 neighbor 外侧
        if (face.neighborCell != -1) {
            const Vector& co = cells_[cellId2index_.at(face.ownerCell)].center;
            const Vector& cn = cells_[cellId2index_.at(face.neighborCell)].center;
            Vector n = face.normal;                   // Face 构造里算出的法向
            Vector d = cn - co;
            if ((n * d) > 0) {                        // 点乘>0 说明朝向进邻居，翻转
                std::vector<int> revIDs = orderedIDs;
                std::reverse(revIDs.begin(), revIDs.end());
                std::reverse(nodeCoords.begin(), nodeCoords.end());
                face = Face(faceId, revIDs, nodeCoords);
                face.ownerCell = kv.second.cells[0];
                face.neighborCell = kv.second.cells.size() == 2 ? kv.second.cells[1] : -1;
            }
        }

        faces_.push_back(face);
        for (int cid : kv.second.cells)
            cells_[cellId2index_.at(cid)].CellFaceIDs.push_back(faceId);
        ++faceId;
    }


	//// === Step 3: 构造 Face 对象并建立 Cell 和 Face 的对应关系 ===
 //   //建立面 点及网格单元之间的拓扑关系
	//unordered_map<vector<int>, vector<int>, VectorHash> FacetoCellMap;  // key:一个面由哪几个节点构成（唯一标识一个面） value:哪些单元拥有这个面（可能是1个或2个）

 //   // 遍历所有单元，提取面信息
 //   for (const auto& cell : cells_) 
 //   {
 //       //遍历组成cell的所有面，并对组成面的Node进行排序，防止重复
 //       for (const auto& faceNodeIDs : cell.getLocalFaces(nodesMap_))  ///注：getLocalFaces()返回的是组成该单元的面，而该面的表示方式是组成该面的两个节点的数组
 //       {
 //           vector<int> sortedNodeIDs = faceNodeIDs;
 //           sort(sortedNodeIDs.begin(), sortedNodeIDs.end()); // 对节点ID排序以确保面的唯一性
	//		FacetoCellMap[sortedNodeIDs].push_back(cell.id); // 将面节点编号作为 key，cell.id 作为 value 存入 faceMap
	//		
 //           //cout << "Cell ID: " << cell.id << " has face with nodes: ";
	//		//for (int nid : sortedNodeIDs) cout << nid << " ";
	//		//cout << endl;
 //       }
 //   }

 //   // 创建面对象并建立与单元的关联
 //   faces_.clear();
 //   int faceId = 1; 
 //   for (const auto& entry : FacetoCellMap)
 //   {//获取网格节点
 //       const vector<int>& nodeIDs = entry.first;  
 //       vector<Vector> nodeCoords;

 //       // 获取面节点的坐标
 //       for (int nid : nodeIDs)
 //           nodeCoords.push_back(nodesMap_[nid].coord); 

 //       Face face(faceId, nodeIDs, nodeCoords); //通过调用构造函数自动计算normal和length
 //       face.ownerCell = entry.second[0];
 //       face.neighborCell = (entry.second.size() == 2) ? entry.second[1] : -1;  //边界网格的边界面只包含Owner网格单元，所有 faceMap的second只有一个元素
 //       
 //       
 //     /*  cout << "构建 Face ID: " << faceId << ", 节点: ";
 //       for (int nid : nodeIDs) cout << nid << " ";*/
 //      // cout << "| ownerCell = " << face.ownerCell
 //          // << ", neighborCell = " << face.neighborCell << endl;

 //       faces_.push_back(face);

 //       for (int cid : entry.second)
 //       {
	//		int idx = cellId2index_.at(cid); //访问 cellId 对应的局部编号 这样才能访问全局编号
 //           cells_[idx].CellFaceIDs.push_back(faceId);
 //           cout << "  → Face ID " << faceId << " 被添加到 Cell ID: " << cid << " 中。" << endl;
 //       }
 //       faceId++;
 //   }
    std::cout << "BuildMesh: faces_ count = " << faces_.size() << std::endl;
    cout << "=== Face owner/neighbor check ===" << std::endl;
    for (const auto& face : faces_) {
       cout << "Face " << face.id << ": [";
        for (int nid : face.FaceNodeIDs)
           cout << nid << " ";
       cout << "], owner = " << face.ownerCell
            << ", neighbor = " << face.neighborCell << std::endl;
    }

	buildFaceBins(); // 构建面 bin
    gmsh::write("UnstructuredMesh-EDFM.msh");
    gmsh::finalize();
}




void Mesh::ClassifySolidMatrixCells()
{
    // 构造一个 faceId 到 Face 的映射
   
    for (auto& cell : cells_)
    {
        bool isBoundary = false;
        for (int faceId : cell.CellFaceIDs)
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

// =========== 计算 (ix,iy) -> 唯一 bin ID  =============
static inline int computeBinID(int ix, int iy, int binCountY)
{
    return ix * binCountY + iy;
}

void Mesh::buildFaceBins()
{
    if (faces_.empty()) return;

    // ---------- 1. 统计整个网格的坐标范围 ----------
    double minX = std::numeric_limits<double>::max();
    double minY = std::numeric_limits<double>::max();
    double maxX = -std::numeric_limits<double>::max();
    double maxY = -std::numeric_limits<double>::max();

    for (const auto& f : faces_) {
        minX = std::min(minX, f.boundingBox.min.m_x);
        minY = std::min(minY, f.boundingBox.min.m_y);
        maxX = std::max(maxX, f.boundingBox.max.m_x);
        maxY = std::max(maxY, f.boundingBox.max.m_y);
    }

    // ---------- 2. 自动计算 bin 数量 ----------
    int targetFacesPerBin = 70;  //可调范围在5-30，越小精度越高，越大过滤速度越快内存消耗越大
    int totalFaces = static_cast<int>(faces_.size());
    int totalBins = std::max(1, totalFaces / targetFacesPerBin);
    int sqrtBins = static_cast<int>(std::sqrt(totalBins));

    binCountX_ = std::min(std::max(sqrtBins, 10), 100); //clamp上下限 5-200 越小 粗网格测试
    binCountY_ = std::min(std::max(sqrtBins, 10), 100);

    binSizeX_ = (maxX - minX) / binCountX_;
    binSizeY_ = (maxY - minY) / binCountY_;

    // ---------- 3. 防护判断 ----------
    if (binSizeX_ <= 0 || binSizeY_ <= 0) {
        std::cerr << "[buildFaceBins] 检测到包围盒尺寸异常，终止构建\n";
        return;
    }

    faceBins_.clear();

    // ---------- 4. 遍历所有面，把它们放进对应 bin ----------
    for (const auto& f : faces_) {
        const AABB& box = f.boundingBox;

        int ix_min = static_cast<int>((box.min.m_x - minX) / binSizeX_);
        int iy_min = static_cast<int>((box.min.m_y - minY) / binSizeY_);
        int ix_max = static_cast<int>((box.max.m_x - minX) / binSizeX_);
        int iy_max = static_cast<int>((box.max.m_y - minY) / binSizeY_);

        ix_min = std::max(0, std::min(binCountX_ - 1, ix_min));
        iy_min = std::max(0, std::min(binCountY_ - 1, iy_min));
        ix_max = std::max(0, std::min(binCountX_ - 1, ix_max));
        iy_max = std::max(0, std::min(binCountY_ - 1, iy_max));

        for (int ix = ix_min; ix <= ix_max; ++ix)
            for (int iy = iy_min; iy <= iy_max; ++iy) {
                int binID = computeBinID(ix, iy, binCountY_);
                faceBins_[binID].push_back(f.id);
            }
    }

    // ---------- 5. 日志 ----------
    std::size_t filledBins = 0;
    for (const auto& kv : faceBins_) if (!kv.second.empty()) ++filledBins;

    std::cout << "[buildFaceBins] 总面数 = " << faces_.size()
        << ", Bin 划分 = " << binCountX_ << "×" << binCountY_
        << ", 平均每 bin 面数 ≈ "
        << static_cast<double>(faces_.size()) / std::max<std::size_t>(1, filledBins)
        << std::endl;
}


/*  ============================================================
 *  Mesh::getCandidateFacesFromBins — 根据查询盒子 (box)
 *  返回可能与之相交的网格面 ID 列表（去重后）。
 *  依赖：binCountX_ / binCountY_ / binSizeX_ / binSizeY_
 *        以及 buildFaceBins() 构建好的 faceBins_。
 *  ============================================================ */
std::vector<int> Mesh::getCandidateFacesFromBins(const AABB& box) const
{
    std::vector<int> result;
    if (faceBins_.empty()) return result;

    // ---------- 1. 全局坐标基准（minX, minY） ----------
    double minX = std::numeric_limits<double>::max();
    double minY = std::numeric_limits<double>::max();
    for (const auto& f : faces_) {
        minX = std::min(minX, f.boundingBox.min.m_x);
        minY = std::min(minY, f.boundingBox.min.m_y);
    }

    // ---------- 2. 计算 AABB 所落入的 bin 区间 ----------
    int ix_min = static_cast<int>((box.min.m_x - minX) / binSizeX_);
    int iy_min = static_cast<int>((box.min.m_y - minY) / binSizeY_);
    int ix_max = static_cast<int>((box.max.m_x - minX) / binSizeX_);
    int iy_max = static_cast<int>((box.max.m_y - minY) / binSizeY_);

    ix_min = std::max(0, std::min(binCountX_ - 1, ix_min));
    iy_min = std::max(0, std::min(binCountY_ - 1, iy_min));
    ix_max = std::max(0, std::min(binCountX_ - 1, ix_max));
    iy_max = std::max(0, std::min(binCountY_ - 1, iy_max));

    // ---------- 3. 预收集候选面 ID（未去重） ----------
    std::vector<int> rawList;
    for (int ix = ix_min; ix <= ix_max; ++ix)
        for (int iy = iy_min; iy <= iy_max; ++iy) {
            int binID = computeBinID(ix, iy, binCountY_);
            auto it = faceBins_.find(binID);
            if (it == faceBins_.end()) continue;
            rawList.insert(rawList.end(), it->second.begin(), it->second.end());
        }

    // ---------- 4. 去重 ----------
    std::sort(rawList.begin(), rawList.end());
    rawList.erase(std::unique(rawList.begin(), rawList.end()), rawList.end());

    // ---------- 5. 二次过滤：AABB overlap ----------
    result.reserve(rawList.size());
    for (int fid : rawList) {
        const Face& f = faces_[fid - 1];  // ID 从 1 开始
        if (f.boundingBox.min.m_x <= box.max.m_x && f.boundingBox.max.m_x >= box.min.m_x &&
            f.boundingBox.min.m_y <= box.max.m_y && f.boundingBox.max.m_y >= box.min.m_y) {
            result.push_back(fid);
        }
    }

    return result;
}

vector<std::reference_wrapper<const Cell>> Mesh::getInnerCells() const
{
    vector<std::reference_wrapper<const Cell>> innerCells;
    for (const auto& cell : cells_)
    {
        if (cell.location == Cell::LocationType::Inner)
            innerCells.push_back(std::cref(cell));
    }
    return innerCells;
}

vector<std::reference_wrapper<const Cell>> Mesh::getBoundaryCells() const
{
	vector<std::reference_wrapper<const Cell>> boundaryCells;
	for (const auto& cell : cells_)
	{
		if (cell.location == Cell::LocationType::Boundary)
			boundaryCells.push_back(std::cref(cell));
	}
	return boundaryCells;
}

void Mesh::printMeshInfo()
{
    std::cout << "总节点数: " << nodes_.size() << endl;
    std::cout << "总单元数: " << cells_.size() << endl;
    std::cout << "总面数: " << faces_.size() << endl;
    int nInner = 0, nBoundary = 0;
    for (const auto& cell : cells_) {
        if (cell.location == Cell::LocationType::Inner) ++nInner;
        else if (cell.location == Cell::LocationType::Boundary) ++nBoundary;
    }
    std::cout << "内部单元数: " << nInner << endl;
    std::cout << "边界单元数: " << nBoundary << endl;

    std::cout << "\n―― 节点信息 ――" << endl;
    for (const auto& node : nodes_)
    {
        std::cout << "Node " << node.id << " : ("
            << node.coord.m_x << ", " << node.coord.m_y << ", " << node.coord.m_z << ")"
            << endl;
    }
    std::cout << "\n―― 单元信息 ――" << endl;      //还需要添加Owner和Neighbor
    for (const auto& cell : cells_)
    {
        if (cell.id < 0) continue;   // 跳过ghost cell
        std::cout << "Cell " << cell.id << " 中心: ("
            << cell.center.m_x << ", " << cell.center.m_y << ", " << cell.center.m_z << ")"
            << " 面积: " << cell.volume << endl;
        std::cout << "  包含面: ";
        for (int fid : cell.CellFaceIDs)
            std::cout << fid << " ";
        std::cout << endl;
    }
    std::cout << "\n―― 面信息 ――" << endl;
    for (const auto& face :faces_)
    {
        std::cout << "Face " << face.id << " 节点: ";
        for (int nid : face.FaceNodeIDs)
            std::cout << nid << " ";
        std::cout << " | owner: " << face.ownerCell
            << " neighbor: " << face.neighborCell << endl;
        std::cout << "  中点: (" << face.midpoint.m_x << ", " << face.midpoint.m_y << ", " << face.midpoint.m_z
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
    for (auto& face : faces_)
    {
        // owner 下标
        auto it_o = cellId2index_.find(face.ownerCell);
        if (it_o == cellId2index_.end()) {
            // 安全：若异常则跳过
            continue;
        }
        const Vector& Cp = cells_[it_o->second].center;

        if (face.neighborCell >= 0)
        {
            // 内部面：用原有 owner+neighbor 的分解
            auto it_n = cellId2index_.find(face.neighborCell);
            if (it_n == cellId2index_.end()) {
                continue;
            }
            const Vector& Cn = cells_[it_n->second].center;
            face.computeFaceVectors(Cp, Cn, method);
        }
        else
        {
            // 边界面：无 ghost，调用新边界分解
            face.computeFaceVectorsBoundary(Cp, method);
        }
    }
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
            << f.FaceNodeIDs[0] << " "
            << f.FaceNodeIDs[1] << " "
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

    //exportGhostInfo(prefix);
}
