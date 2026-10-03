#pragma once
#include <cmath>
#include <stdexcept>

// Physical channel dimensions in meters, immutable for the lifetime of a run.
struct channel_dimensions {
    const double W, H, L;
    explicit channel_dimensions(double width = 5e-4, double height = 4e-5, double length = 0.025)
        : W(width), H(height), L(length) {
        if (!std::isfinite(W) || !std::isfinite(H) || !std::isfinite(L) || W <= 0 || H <= 0 || L <= 0)
            throw std::invalid_argument("Channel dimensions must be finite positive values in meters.");
    }
};
