/**
 * @file ADVar.hpp
 * @brief 轻量级前向自动微分 (Forward-Mode Automatic Differentiation) 核心库
 * @details 专为全隐式方法 (FIM) 设计，基于 Eigen 库存储局部梯度。
 * 实现了所有基础算数运算符和超越函数的重载，自动且精确地应用链式法则。
 */

#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <limits>

 /**
  * @class ADVar
  * @brief 前向自动微分变量类
  * @tparam N 独立自变量的数量（即每个网格块的自由度 DOF，如单相 N=2, 两相 N=3）
  */
template<int N>
class ADVar
{
public:
    double val;                               ///< 变量的物理值
    Eigen::Matrix<double, N, 1> grad;         ///< 变量对主元 (P, S, T 等) 的偏导数向量

    // =========================================================
    // 1. 构造函数与初始化
    // =========================================================

    /**
     * @brief 默认构造函数，初始化为 0.0，梯度为零向量
     */
    ADVar() : val(0.0)
    {
        grad.setZero();
    }

    /**
     * @brief 常数构造函数
     * @param v 常数值 (其梯度自动设为 0)
     */
    ADVar(double v) : val(v)
    {
        grad.setZero();
    }

    /**
     * @brief 独立自变量构造函数 (主变量初始化)
     * @details 用于将该变量标记为求导基准元。
     * 例如：初始化压力 P，使其对自身导数为 1 (即 grad[idx] = 1.0)
     * @param v 变量初始值
     * @param idx 自变量在系统中的局部自由度索引 (0 <= idx < N)
     */
    ADVar(double v, int idx) : val(v)
    {
        grad.setZero();
        if (idx >= 0 && idx < N)
        {
            grad(idx) = 1.0;
        }
        else
        {
            std::cerr << "[ADVar Error] Initialization index " << idx << " out of bounds [0, " << N - 1 << "]." << std::endl;
        }
    }

    // =========================================================
    // 2. 状态检查与安全防护 (Safety & Validation)
    // =========================================================

    /**
     * @brief 检查当前变量是否包含 NaN 或 Inf
     * @return bool 如果包含无效数值则返回 true
     */
    bool hasNaN() const
    {
        if (std::isnan(val) || std::isinf(val)) return true;
        for (int i = 0; i < N; ++i)
        {
            if (std::isnan(grad(i)) || std::isinf(grad(i))) return true;
        }
        return false;
    }

    // =========================================================
    // 3. 复合赋值运算符重载 (Compound Assignment)
    // =========================================================

    ADVar<N>& operator+=(const ADVar<N>& rhs)
    {
        val += rhs.val;
        grad += rhs.grad;
        return *this;
    }

    ADVar<N>& operator+=(double rhs)
    {
        val += rhs;
        return *this;
    }

    ADVar<N>& operator-=(const ADVar<N>& rhs)
    {
        val -= rhs.val;
        grad -= rhs.grad;
        return *this;
    }

    ADVar<N>& operator-=(double rhs)
    {
        val -= rhs;
        return *this;
    }

    ADVar<N>& operator*=(const ADVar<N>& rhs)
    {
        // 乘法链式法则: d(u*v) = v*du + u*dv
        grad = rhs.val * grad + val * rhs.grad;
        val *= rhs.val;
        return *this;
    }

    ADVar<N>& operator*=(double rhs)
    {
        val *= rhs;
        grad *= rhs;
        return *this;
    }

    ADVar<N>& operator/=(const ADVar<N>& rhs)
    {
        // 除法链式法则: d(u/v) = (v*du - u*dv) / v^2
        if (std::abs(rhs.val) < 1e-30)
        {
            std::cerr << "[ADVar Critical] Division by near-zero value during /= operation!" << std::endl;
        }
        double inv_v = 1.0 / rhs.val;
        double inv_v2 = inv_v * inv_v;
        grad = (grad * rhs.val - val * rhs.grad) * inv_v2;
        val *= inv_v;
        return *this;
    }

    ADVar<N>& operator/=(double rhs)
    {
        if (std::abs(rhs) < 1e-30)
        {
            std::cerr << "[ADVar Critical] Division by near-zero scalar during /= operation!" << std::endl;
        }
        double inv_rhs = 1.0 / rhs;
        val *= inv_rhs;
        grad *= inv_rhs;
        return *this;
    }

    // 一元负号
    ADVar<N> operator-() const
    {
        ADVar<N> res;
        res.val = -val;
        res.grad = -grad;
        return res;
    }
};

// =========================================================
// 4. 二元算术运算符重载 (Binary Arithmetic Operators)
// =========================================================

template<int N> inline ADVar<N> operator+(const ADVar<N>& lhs, const ADVar<N>& rhs) {
    ADVar<N> res = lhs; return res += rhs;
}
template<int N> inline ADVar<N> operator+(const ADVar<N>& lhs, double rhs) {
    ADVar<N> res = lhs; return res += rhs;
}
template<int N> inline ADVar<N> operator+(double lhs, const ADVar<N>& rhs) {
    ADVar<N> res = rhs; return res += lhs; // 加法交换律
}

template<int N> inline ADVar<N> operator-(const ADVar<N>& lhs, const ADVar<N>& rhs) {
    ADVar<N> res = lhs; return res -= rhs;
}
template<int N> inline ADVar<N> operator-(const ADVar<N>& lhs, double rhs) {
    ADVar<N> res = lhs; return res -= rhs;
}
template<int N> inline ADVar<N> operator-(double lhs, const ADVar<N>& rhs) {
    ADVar<N> res;
    res.val = lhs - rhs.val;
    res.grad = -rhs.grad;
    return res;
}

template<int N> inline ADVar<N> operator*(const ADVar<N>& lhs, const ADVar<N>& rhs) {
    ADVar<N> res = lhs; return res *= rhs;
}
template<int N> inline ADVar<N> operator*(const ADVar<N>& lhs, double rhs) {
    ADVar<N> res = lhs; return res *= rhs;
}
template<int N> inline ADVar<N> operator*(double lhs, const ADVar<N>& rhs) {
    ADVar<N> res = rhs; return res *= lhs; // 乘法交换律
}

template<int N> inline ADVar<N> operator/(const ADVar<N>& lhs, const ADVar<N>& rhs) {
    ADVar<N> res = lhs; return res /= rhs;
}
template<int N> inline ADVar<N> operator/(const ADVar<N>& lhs, double rhs) {
    ADVar<N> res = lhs; return res /= rhs;
}
template<int N> inline ADVar<N> operator/(double lhs, const ADVar<N>& rhs) {
    if (std::abs(rhs.val) < 1e-30) {
        std::cerr << "[ADVar Critical] Division by near-zero value in (double / ADVar)!" << std::endl;
    }
    ADVar<N> res;
    double inv_v = 1.0 / rhs.val;
    res.val = lhs * inv_v;
    res.grad = (-lhs * inv_v * inv_v) * rhs.grad;
    return res;
}

// =========================================================
// 5. 核心数学超越函数重载 (Transcendental Math Functions)
// =========================================================

/**
 * @brief 指数函数 exp(u)
 * @details d(exp(u)) = exp(u) * du
 */
template<int N> inline ADVar<N> exp(const ADVar<N>& u)
{
    ADVar<N> res;
    res.val = std::exp(u.val);
    res.grad = res.val * u.grad;
    return res;
}

/**
 * @brief 自然对数函数 log(u)
 * @details d(log(u)) = (1/u) * du
 */
template<int N> inline ADVar<N> log(const ADVar<N>& u)
{
    if (u.val <= 0.0) {
        std::cerr << "[ADVar Critical] Logarithm of non-positive value: " << u.val << std::endl;
    }
    ADVar<N> res;
    res.val = std::log(u.val);
    res.grad = (1.0 / u.val) * u.grad;
    return res;
}

/**
 * @brief 幂函数 pow(u, a)
 * @details d(u^a) = a * u^(a-1) * du
 */
template<int N> inline ADVar<N> pow(const ADVar<N>& u, double a)
{
    ADVar<N> res;
    res.val = std::pow(u.val, a);
    if (u.val == 0.0 && a < 1.0) {
        std::cerr << "[ADVar Critical] Derivative singularity in pow(0, a) where a < 1." << std::endl;
    }
    double deriv = a * std::pow(u.val, a - 1.0);
    res.grad = deriv * u.grad;
    return res;
}

/**
 * @brief 平方根函数 sqrt(u)
 * @details d(sqrt(u)) = 0.5 / sqrt(u) * du
 */
template<int N> inline ADVar<N> sqrt(const ADVar<N>& u)
{
    if (u.val < 0.0) {
        std::cerr << "[ADVar Critical] Square root of negative value: " << u.val << std::endl;
    }
    ADVar<N> res;
    res.val = std::sqrt(u.val);

    // 防范 x=0 时导数为无穷大的奇点
    double deriv = (res.val > 1e-30) ? (0.5 / res.val) : 0.0;
    res.grad = deriv * u.grad;
    return res;
}

/**
 * @brief 绝对值函数 abs(u)
 * @details d(|u|) = sign(u) * du. (在 u=0 处定义导数为 0 以保数值稳定)
 */
template<int N> inline ADVar<N> abs(const ADVar<N>& u)
{
    ADVar<N> res;
    res.val = std::abs(u.val);
    double sign = (u.val > 0.0) ? 1.0 : ((u.val < 0.0) ? -1.0 : 0.0);
    res.grad = sign * u.grad;
    return res;
}

/**
 * @brief 正弦函数 sin(u)
 */
template<int N> inline ADVar<N> sin(const ADVar<N>& u)
{
    ADVar<N> res;
    res.val = std::sin(u.val);
    res.grad = std::cos(u.val) * u.grad;
    return res;
}

/**
 * @brief 余弦函数 cos(u)
 */
template<int N> inline ADVar<N> cos(const ADVar<N>& u)
{
    ADVar<N> res;
    res.val = std::cos(u.val);
    res.grad = -std::sin(u.val) * u.grad;
    return res;
}

// =========================================================
// 6. 比较运算符重载 (Logical Operators)
// =========================================================
// 注：逻辑比较仅比较物理值 (val)，用于条件分支判断 (如迎风格式的流向判断)

template<int N> inline bool operator<(const ADVar<N>& lhs, const ADVar<N>& rhs) { return lhs.val < rhs.val; }
template<int N> inline bool operator<(const ADVar<N>& lhs, double rhs) { return lhs.val < rhs; }
template<int N> inline bool operator<(double lhs, const ADVar<N>& rhs) { return lhs < rhs.val; }

template<int N> inline bool operator>(const ADVar<N>& lhs, const ADVar<N>& rhs) { return lhs.val > rhs.val; }
template<int N> inline bool operator>(const ADVar<N>& lhs, double rhs) { return lhs.val > rhs; }
template<int N> inline bool operator>(double lhs, const ADVar<N>& rhs) { return lhs > rhs.val; }

template<int N> inline bool operator<=(const ADVar<N>& lhs, const ADVar<N>& rhs) { return lhs.val <= rhs.val; }
template<int N> inline bool operator<=(const ADVar<N>& lhs, double rhs) { return lhs.val <= rhs; }
template<int N> inline bool operator<=(double lhs, const ADVar<N>& rhs) { return lhs <= rhs.val; }

template<int N> inline bool operator>=(const ADVar<N>& lhs, const ADVar<N>& rhs) { return lhs.val >= rhs.val; }
template<int N> inline bool operator>=(const ADVar<N>& lhs, double rhs) { return lhs.val >= rhs; }
template<int N> inline bool operator>=(double lhs, const ADVar<N>& rhs) { return lhs >= rhs.val; }

// =========================================================
// 7. Max / Min 辅助函数 (非平滑函数的局部次导数)
// =========================================================

template<int N> inline ADVar<N> max(const ADVar<N>& lhs, const ADVar<N>& rhs) {
    return (lhs.val > rhs.val) ? lhs : rhs;
}
template<int N> inline ADVar<N> max(const ADVar<N>& lhs, double rhs) {
    return (lhs.val > rhs) ? lhs : ADVar<N>(rhs); // 常数梯度为 0
}
template<int N> inline ADVar<N> max(double lhs, const ADVar<N>& rhs) {
    return (lhs > rhs.val) ? ADVar<N>(lhs) : rhs;
}

template<int N> inline ADVar<N> min(const ADVar<N>& lhs, const ADVar<N>& rhs) {
    return (lhs.val < rhs.val) ? lhs : rhs;
}
template<int N> inline ADVar<N> min(const ADVar<N>& lhs, double rhs) {
    return (lhs.val < rhs) ? lhs : ADVar<N>(rhs);
}
template<int N> inline ADVar<N> min(double lhs, const ADVar<N>& rhs) {
    return (lhs < rhs.val) ? ADVar<N>(lhs) : rhs;
}