pub const Backend = enum { chebyshev, newton };

pub const Interval = struct {
    t0: f64,
    dt: f64,
    // convenience:
    pub fn toY(self: Interval, t: f64) f64 {
        // map t -> y in [-1,1]
        _ = .{ self, t };
    }
    pub fn toT(self: Interval, y: f64) f64 {
        // map y -> t
        _ = .{ self, y };
    }
};

pub const ChebyshevGrid = struct {
    // boundary-including nodes y_l = -cos(l*pi/(M-1))  (Eq. 76, J. Comp. Phys. 343, 368, 2017)
    pub fn nodesBoundaryIncluding(alloc: *mem.Allocator, M: usize) ![]f64 {
        _ = .{ alloc, M };
    }
    pub fn mapToInterval(alloc: *mem.Allocator, y: []const f64, interval: Interval) ![]f64 {
        _ = .{ alloc, y, interval };
    }
};

pub fn Polynomial(comptime T: type, comptime backend: Backend) type {
    _ = .{ T, backend };

    return struct {
        pub const Coeffs = struct {
            // for chebyshev: c[0..M-1]
            // for newton:    a[0..M-1]
            data: []T,
            // explicit: only used by newton-stability transform (length-4 mapping)
            domain_scale: f64 = 1.0,
        };

        pub fn build(alloc: *mem.Allocator, y_nodes: []const f64, f: []const T, interval: Interval) !Coeffs {
            _ = .{ alloc, y_nodes, f, interval };
        }

        // evaluate at normalized y (can be slightly outside [-1,1] if caller wants extrapolation)
        pub fn eval(y: f64, coeffs: Coeffs, y_nodes: []const f64) T {
            _ = .{ y, coeffs, y_nodes };
        }
    };
}

const std = @import("std");
const mem = std.mem;
