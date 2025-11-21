#pragma once

namespace IMPES
{
    template<typename T>
    inline T clampValue(T v, T lo, T hi)
    {
        if (v < lo) return lo;
        if (v > hi) return hi;
        return v;
    }
} // namespace IMPES

namespace IMPES_revised
{
    template<typename T>
    inline T clampValue(T v, T lo, T hi)
    {
        if (v < lo) return lo;
        if (v > hi) return hi;
        return v;
    }

}

