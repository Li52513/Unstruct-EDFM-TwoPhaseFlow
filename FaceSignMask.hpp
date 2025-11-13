#pragma once

#include <cstdint>
#include <vector>
#include <cmath>

#include "FaceFieldRegistry.h"

struct FaceSignMask
{
    std::vector<int8_t> signBits;
    bool initialized = false;

    void ensureSize(size_t n)
    {
        if (signBits.size() != n) {
            signBits.assign(n, 0);
            initialized = false;
        }
    }

    void reset()
    {
        signBits.clear();
        initialized = false;
    }
};

struct FaceSignUpdateInfo
{
    bool hadPreviousMask = false;
    size_t flippedCount = 0;
    size_t totalFaces = 0;

    double flipRatio() const
    {
        return (totalFaces == 0) ? 0.0 : static_cast<double>(flippedCount) / static_cast<double>(totalFaces);
    }
};

inline FaceSignUpdateInfo updateFaceSignMask_fromFlux(
    const faceScalarField& flux,
    double epsilon,
    FaceSignMask& mask,
    std::vector<int>& flippedFaces)
{
    FaceSignUpdateInfo info;
    const size_t nFaces = flux.data.size();
    info.totalFaces = nFaces;

    mask.ensureSize(nFaces);
    flippedFaces.clear();
    info.hadPreviousMask = mask.initialized;

    for (size_t i = 0; i < nFaces; ++i) {
        const double val = flux.data[i];
        int8_t newSign = 0;
        if (std::abs(val) > epsilon) {
            newSign = (val > 0.0) ? 1 : -1;
        }

        if (!mask.initialized) {
            mask.signBits[i] = newSign;
            continue;
        }

        if (newSign != mask.signBits[i]) {
            mask.signBits[i] = newSign;
            flippedFaces.push_back(static_cast<int>(i));
        }
    }

    if (!mask.initialized) {
        mask.initialized = true;
    }
    info.flippedCount = flippedFaces.size();
    return info;
}

