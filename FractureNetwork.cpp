#include "FractureNetwork.h"
#include "8_DOP.h"
#include <algorithm>
#include <cmath>
#include <iostream> 
#include <iomanip>  
#include <unordered_map>
#include <tuple>
#include <numeric>      // accumulate
#include <cmath>

// 自定义 PI
static constexpr double PI = 3.14159265358979323846;

///------------DFN生成相关辅助函数-------------///
// uniform [0,1)
static double uniformRand() {
    return double(std::rand()) / double(RAND_MAX);
}

// von Mises (Best–Fisher), mean=0
static double sampleVonMises(double kappa) {
    if (kappa < 1e-6) {
        return uniformRand() * 2.0 * PI;
    }
    double a = 1.0 + std::sqrt(1.0 + 4.0 * kappa * kappa);
    double b = (a - std::sqrt(2.0 * a)) / (2.0 * kappa);
    double r = (1.0 + b * b) / (2.0 * b);
    while (true) {
        double u1 = uniformRand();
        double u2 = uniformRand();
        double z = std::cos(PI * u1);
        double f = (1.0 + r * z) / (r + z);
        double c = kappa * (r - f);
        if (u2 < c * (2.0 - c) || std::log(u2) <= c - 1.0) 
        {
            double u3 = uniformRand();
            double theta = std::acos(f);
            return (u3 > 0.5 ? theta : 2.0 * PI - theta);
        }
    }
}

// 2D segment–segment intersection test
static bool segmentsIntersect(const Vector& a1, const Vector& a2,
    const Vector& b1, const Vector& b2) {
    auto cross = [](double x1, double y1, double x2, double y2) {
        return x1 * y2 - y1 * x2;
        };
    Vector r{ a2.m_x - a1.m_x, a2.m_y - a1.m_y, 0.0 };
    Vector s{ b2.m_x - b1.m_x, b2.m_y - b1.m_y, 0.0 };
    double denom = cross(r.m_x, r.m_y, s.m_x, s.m_y);
    if (std::fabs(denom) < 1e-12) return false;
    Vector d{ b1.m_x - a1.m_x, b1.m_y - a1.m_y, 0.0 };
    double t = cross(d.m_x, d.m_y, s.m_x, s.m_y) / denom;
    double u = cross(d.m_x, d.m_y, r.m_x, r.m_y) / denom;
    return (t >= 0.0 && t <= 1.0 && u >= 0.0 && u <= 1.0);
}

void FractureNetwork::setRandomSeed(unsigned seed)
{
	std::srand(seed);
}
// ----------------------------------------------------------------------------


void FractureNetwork::generateDFN(
    int N,
    const Vector& minPoint,
    const Vector& maxPoint,
    double Lmin,
    double Lmax,
    double alpha,
    double kappa,
    bool avoidOverlap
) {
    // helper: power-law length
    auto sampleLength = [&]() {
        double u = uniformRand();
        double c0 = std::pow(Lmin, 1.0 - alpha);
        double c1 = std::pow(Lmax, 1.0 - alpha);
        return std::pow(c0 + u * (c1 - c0), 1.0 / (1.0 - alpha));
        };

    int created = 0, attempts = 0, maxAttempts = N * 10;
    while (created < N && attempts < maxAttempts) {
        ++attempts;
        // 1) center
        double cx = minPoint.m_x + uniformRand() * (maxPoint.m_x - minPoint.m_x);
        double cy = minPoint.m_y + uniformRand() * (maxPoint.m_y - minPoint.m_y);
        double cz = minPoint.m_z + uniformRand() * (maxPoint.m_z - minPoint.m_z);
        Vector C{ cx,cy,cz };

        // 2) length & orientation
        double L = sampleLength();
        double theta = sampleVonMises(kappa);

        // 3) endpoints
        double dx = 0.5 * L * std::cos(theta),
            dy = 0.5 * L * std::sin(theta);
        Vector s{ C.m_x + dx, C.m_y + dy, C.m_z };
        Vector e{ C.m_x - dx, C.m_y - dy, C.m_z };

        // 4) boundary check
        if (s.m_x<minPoint.m_x || s.m_x>maxPoint.m_x ||
            s.m_y<minPoint.m_y || s.m_y>maxPoint.m_y ||
            e.m_x<minPoint.m_x || e.m_x>maxPoint.m_x ||
            e.m_y<minPoint.m_y || e.m_y>maxPoint.m_y)
            continue;

        // 5) optional overlap avoidance
        if (avoidOverlap) {
            bool bad = false;
            for (auto const& F : fractures) {
                if (!F.intersections.empty()) {
                    Vector a1 = F.intersections.front().point;
                    Vector a2 = F.intersections.back().point;
                    if (segmentsIntersect(s, e, a1, a2)) { bad = true; break; }
                }
            }
            if (bad) continue;
        }

        // 6) finally add
        addFracture(s, e);
        ++created;
    }

    if (created < N) {
        std::cerr << "[Warning] only " << created
            << " of " << N << " fractures.\n";
    }
}

//-------------------------------------------------------------------------------------------------//

//-------------添加一条裂缝--------------//

void FractureNetwork::addFracture(const Vector& s, const Vector& e)
{
    fractures.emplace_back(s, e);           // 1. 先放进容器
    auto& f = fractures.back();             // 2. 取出刚刚插入的那条
    f.id = nextFracID_++;                // 3. 写入唯一 ID
}

/// 检测裂缝之间的交点
void FractureNetwork::DetectFracturetoFractureIntersections()
{
    globalFFPts.clear();
    int nextID = 1;
    for (size_t i = 0; i < fractures.size(); ++i) 
    {
        const auto& Fi = fractures[i];
        Vector p = Fi.start, q = Fi.end;
        double L = (q - p).Mag();
        //cout << "第" << Fi.id << "条裂缝的长度为" << L << endl;
        for (size_t j = i + 1; j < fractures.size(); ++j) 
        {
            const auto& Fj = fractures[j];
            Vector r = Fj.start, s = Fj.end;
            Vector ip;
            if (Fracture::lineSegmentIntersection(p, q, r, s, ip)) 
            {
                //cout << "第" << Fi.id << "条和第" << Fj.id << "条裂缝的交点坐标为（" << ip.m_x << "," << ip.m_y << "," << ip.m_z << "）" << endl;
                // 计算各自的 param
                double ti = L > 1e-12 ? (ip - p).Mag() / L : 0.0;
                double L2 = (s - r).Mag();
                double tj = L2 > 1e-12 ? (ip - r).Mag() / L2 : 0.0;
                globalFFPts.push_back({ nextID++, ip, int(i), int(j), ti, tj });
            }
        }
    }
}

/// 对裂缝与裂缝的交点进行去重并编号
void FractureNetwork::DeduplicateAndRenumberFractureToFractureIntersections()
{
    vector<GlobalFFPoint> uniquePts; //
    int gid = 1;
    for (auto& g : globalFFPts)
    {
        bool found = false;
        for (auto& u : uniquePts)
        {
            if (isClose(g.point, u.point)) { found = true; break; }
        }
        if (!found)
        {
            g.id = gid++;
            uniquePts.push_back(g);
        }
    }
    globalFFPts.swap(uniquePts);
}


void FractureNetwork::DistributeFracture_FractureIntersectionsToGlobalInersections()
{
    constexpr double tol = 1e-8;
    for (auto& g : globalFFPts)
    {
        // —— 裂缝 A —— 
        {
            auto& FA = fractures[g.fracA];
            // 在 FA.intersections 中查找：有没有已有的交点其 param 与 g.paramA 极其接近
            auto it = std::find_if(FA.intersections.begin(), FA.intersections.end(),
                [&](auto& ip) {
                    return std::abs(ip.param - g.paramA) < tol;
                });
            if (it != FA.intersections.end()) {
                // 找到了，就“升级”这条记录，把它标记为裂缝–裂缝交点
                it->isFF = true;
                it->globalFFID = g.id;
                it->origin = IntersectionOrigin::FracFrac;
            }
            else {
                // 没找到，就新 emplace 一个节点
                FA.intersections.emplace_back(
                    -1,                   // id 先用 -1，后面统一重编号
                    g.point,              // 坐标：交点位置
                    -1,                   // edgeID：不是与网格面相交
                    g.paramA,             // param：在裂缝 A 上的位置归一化参数
                    true,                 // isFF：确实是裂缝–裂缝交点
                    g.id,                 // globalFFID：全局交点编号
                    IntersectionOrigin::FracFrac  // origin：来自“裂缝–裂缝”
                );
            }
        }

        // —— 裂缝 B —— 完全同理，只不过用 g.paramB，fractures[g.fracB]
        {
            auto& FB = fractures[g.fracB];
            auto it = std::find_if(FB.intersections.begin(), FB.intersections.end(),
                [&](auto& ip) {
                    return std::abs(ip.param - g.paramB) < tol;
                });
            if (it != FB.intersections.end()) {
                it->isFF = true;
                it->globalFFID = g.id;
                it->origin = IntersectionOrigin::FracFrac;
            }
            else {
                FB.intersections.emplace_back(
                    -1,
                    g.point,
                    -1,
                    g.paramB,
                    true,
                    g.id,
                    IntersectionOrigin::FracFrac
                );
            }
        }
    }
}


/// 判断两个点是否在给定的公差范围内相等
bool FractureNetwork::isClose(const Vector& a, const Vector& b, double tol) const
{
    return (a - b).Mag() < tol;
}


///
const FractureElement* FractureNetwork::getElementByGlobalID(int globalID) const
{
    if (!fracElemIndexValid_) return nullptr;
    const auto& offsets = fracElemIndex_.offset;

    // 1. 越界检查
    if (globalID < 0 || globalID >= (int)fracElemIndex_.total) return nullptr;

    // 2. 二分查找：找到 globalID 落在哪个 offset 区间
    // offsets 是递增序列，例如 [0, 5, 8, 12 ...]
    // upper_bound 返回第一个 > globalID 的位置
    auto it = std::upper_bound(offsets.begin(), offsets.end(), globalID);

    // 对应的 Fracture 索引是 (it - begin) - 1
    int fracID = std::distance(offsets.begin(), it) - 1;

    // 3. 计算局部索引
    // Local 0-based index = GlobalID - Offset[fracID]
    int localIndex = globalID - offsets[fracID];

    // 4. 返回对象
    // fractures[fracID].elements 是 vector，0-based，直接对应 localIndex
    if (fracID >= 0 && fracID < (int)fractures.size())
    {
        const auto& elems = fractures[fracID].elements;
        if (localIndex >= 0 && localIndex < (int)elems.size())
        {
            return &elems[localIndex];
        }
    }
    return nullptr;

}


void FractureNetwork::printFractureInfo() const
{
    std::cout << "\n========= Fracture Information =========\n";
    for (size_t fid = 0; fid < fractures.size(); ++fid)
    {
        const auto& F = fractures[fid];
        std::cout << "Fracture #" << fid + 1 << "  ("
            << F.start.m_x << "," << F.start.m_y << ")  →  ("
            << F.end.m_x << "," << F.end.m_y << ")\n";
    }

    /* 全局 FF 交点一览 ---------------------------------------------------*/
    std::cout << "\n========= Global Fracture–Fracture Intersections =========\n";
    for (const auto& g : globalFFPts)
    {
        std::cout << "ID " << g.id << "  :  Frac "
            << g.fracA + 1 << "  &  " << g.fracB + 1
            << "   (" << g.point.m_x << "," << g.point.m_y << ")\n";
    }
}

void FractureNetwork::exportToTxt(const std::string& prefix) const
{
    /*================ 1. 裂缝段端点文件  ================*/
    std::ofstream fsSeg(prefix + "_fractureSegments.txt");
    fsSeg << "fid segid  x1  y1  x2  y2  mx  my  cellID\n";

    for (size_t fid = 0; fid < fractures.size(); ++fid)
    {
        const auto& F = fractures[fid];

        for (const auto& elem : F.elements)
        {
            /* 找到该段在 intersections 中的两端点：
               subdivide() 保证 elem.id == its position starting at 1 */
            const auto& I = F.intersections;
            const Vector& p1 = I[elem.id - 1].point;
            const Vector& p2 = I[elem.id].point;

            Vector mid = 0.5 * (p1 + p2);

            fsSeg << fid << " "
                << elem.id << " "
                << p1.m_x << " " << p1.m_y << " "
                << p2.m_x << " " << p2.m_y << " "
                << mid.m_x << " " << mid.m_y << " "
                << elem.cellID << "\n";
        }
    }
    fsSeg.close();

    /*================ 2. 裂缝-裂缝交点坐标  ================*/
    std::ofstream fsFF(prefix + "_ffPoints.txt");
    fsFF << "id  x  y  fracA  segA  fracB  segB\n";

    /* 为避免星形交点重复，两点容差以内只写一次。
       这里直接用 globalFFPts 的 id，若同一点有多对 fracA/fracB，
       会生成多条记录，但坐标相同；可由绘图脚本自动合并显示。   */

    for (const auto& g : globalFFPts)
    {
        fsFF << g.id << " "
            << g.point.m_x << " " << g.point.m_y << " "
            << g.fracA + 1 << " "
            << fractures[g.fracA].locateSegment(g.paramA) + 1 << " "
            << g.fracB + 1 << " "
            << fractures[g.fracB].locateSegment(g.paramB) + 1 << "\n";
    }
    fsFF.close();

    /*================ 3. 可选：TI 系数  ===================*/
    std::ofstream fsTI(prefix + "_ffTI.csv");
    fsTI << "fracA,segA,fracB,segB,TIw,TIg\n";
    for (size_t fid = 0; fid < fractures.size(); ++fid) {
        const auto& F = fractures[fid];
        for (size_t sid = 0; sid < F.elements.size(); ++sid) {
            const auto& elem = F.elements[sid];
            for (const auto& ex : elem.ffEx) {
                // 只写一次（避免重复）：
                if (fid < (size_t)ex.peerFracID || (fid == (size_t)ex.peerFracID && sid < (size_t)ex.peerSegID)) {
                    fsTI << (fid + 1) << "," << (sid + 1) << ","
                        << (ex.peerFracID + 1) << "," << (ex.peerSegID + 1) << ","
                        << ex.TIw << "," << ex.TIg << "\n";
                }
            }
        }
    }
    fsTI.close();

    /*================ 4. [新增] 导出 8-DOP 几何尺寸数据 ================*/
    // 格式: fid minX maxX minY maxY minD1 maxD1 minD2 maxD2
    std::string filename8DOP = prefix + "_8DOP.txt";
    std::ofstream fs8DOP(filename8DOP);

    if (!fs8DOP.is_open()) {
        std::cerr << "[Error] Unable to create 8-DOP export file: " << filename8DOP << std::endl;
    }
    else {
        // 输出表头
        fs8DOP << "fid minX maxX minY maxY minD1 maxD1 minD2 maxD2\n";

        for (const auto& F : fractures)
        {
            // 实时构建 8-DOP (复用已验证的逻辑)
            Box8DOP box;
            box.fromSegment(F.start, F.end);

            fs8DOP << F.id << " "
                << box.minX << " " << box.maxX << " "
                << box.minY << " " << box.maxY << " "
                << box.minD1 << " " << box.maxD1 << " " // D1 = x + y
                << box.minD2 << " " << box.maxD2 << "\n"; // D2 = x - y
        }
        fs8DOP.close();
        std::cout << "[Export] 8-DOP data exported to " << filename8DOP << std::endl;
    }

    /*================ 5. [新增] 导出 AABB 几何尺寸数据 ================*/
    // 格式: fid minX maxX minY maxY
    std::string filenameAABB = prefix + "_AABB.txt";
    std::ofstream fsAABB(filenameAABB);

    if (!fsAABB.is_open()) {
        std::cerr << "[Error] Unable to create AABB export file: " << filenameAABB << std::endl;
    }
    else {
        // 输出表头
        fsAABB << "fid minX maxX minY maxY\n";

        for (const auto& F : fractures)
        {
            // 直接利用起点终点计算，确保数据绝对准确（避免依赖可能未更新的缓存状态）
            double minX = std::min(F.start.m_x, F.end.m_x);
            double maxX = std::max(F.start.m_x, F.end.m_x);
            double minY = std::min(F.start.m_y, F.end.m_y);
            double maxY = std::max(F.start.m_y, F.end.m_y);

            fsAABB << F.id << " "
                << minX << " " << maxX << " "
                << minY << " " << maxY << "\n";
        }
        fsAABB.close();
        std::cout << "[Export] AABB data exported to " << filenameAABB << std::endl;
    }
}