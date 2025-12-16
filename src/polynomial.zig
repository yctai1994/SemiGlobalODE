const Backend = enum { chebyshev, newton };

const Interval = struct {
    t0: f64,
    dt: f64,

    const IntervalError = error{InvalidInterval};

    fn init(t0: f64, dt: f64) !Interval {
        if (dt <= 0.0) return IntervalError.InvalidInterval;
        return .{ .t0 = t0, .dt = dt };
    }

    // map t -> y in [-1,1]
    inline fn spaceToChebyshev(self: *const Interval, t: f64) f64 {
        return (2.0 / self.dt) * (t - self.t0) - 1;
    }

    // map y -> t
    inline fn chebyshevToSpace(self: *const Interval, y: f64) f64 {
        return self.t0 + 0.5 * self.dt * (1 + y);
    }
};

const GridError = error{
    InvalidOrder,
    LengthMismatch,
};

/// "Chebyshev" means the dimensionless reference coordinate y ∈ [-1,1]
const ChebyshevGrid = struct {
    grid: []f64,

    // boundary-including nodes y(l) = -cos( l*π / (M-1) )  (Eq. 76, J. Comp. Phys. 343, 368, 2017)
    fn init(allocator: mem.Allocator, order: usize) !*ChebyshevGrid {
        if (order < 2) return GridError.InvalidOrder;

        const self: *ChebyshevGrid = try allocator.create(ChebyshevGrid);
        errdefer allocator.destroy(self);

        self.grid = try allocator.alloc(f64, order);
        errdefer allocator.free(self.grid);

        const de: f64 = @floatFromInt(order - 1);
        for (self.grid, 0..) |*ptr, ind| {
            const nu: f64 = @floatFromInt(ind);
            ptr.* = -@cos(math.pi * (nu / de));
        }

        return self;
    }

    fn deinit(self: *ChebyshevGrid, allocator: mem.Allocator) void {
        allocator.free(self.grid);
        allocator.destroy(self);
    }

    fn chebyshevToSpace(chebyshev: *const ChebyshevGrid, space: *SpaceGrid, interval: Interval) !void {
        if (chebyshev.grid.len != space.grid.len) return GridError.LengthMismatch;
        for (space.grid, chebyshev.grid) |*dest, src| dest.* = interval.chebyshevToSpace(src);
    }
};

/// "Space" means the physical independent variable (t, x, etc.)
const SpaceGrid = struct {
    grid: []f64,

    fn initUninitialized(allocator: mem.Allocator, order: usize) !*SpaceGrid {
        if (order < 2) return GridError.InvalidOrder;

        const self: *SpaceGrid = try allocator.create(SpaceGrid);
        errdefer allocator.destroy(self);

        self.grid = try allocator.alloc(f64, order);
        errdefer allocator.free(self.grid);

        return self;
    }

    fn deinit(self: *SpaceGrid, allocator: mem.Allocator) void {
        allocator.free(self.grid);
        allocator.destroy(self);
    }

    fn spaceToChebyshev(space: *const SpaceGrid, chebyshev: *ChebyshevGrid, interval: Interval) !void {
        if (space.grid.len != chebyshev.grid.len) return GridError.LengthMismatch;
        for (chebyshev.grid, space.grid) |*dest, src| dest.* = interval.spaceToChebyshev(src);
    }
};

/// We'll likely have lots of sampled things (u( t[l] ), G( t[l] ), s( t[l] ), etc.)
const GridSamples = struct {
    // owns: f( t[l] ) values
    vals: []f64,

    // t[l] values, does not own;
    // just references to a SpaceGrid which must outlive GridSamples
    grid: *const SpaceGrid,

    fn init(allocator: mem.Allocator, space: *const SpaceGrid) !*GridSamples {
        const self: *GridSamples = try allocator.create(GridSamples);
        errdefer allocator.destroy(self);

        self.vals = try allocator.alloc(f64, space.grid.len);
        errdefer allocator.free(self.vals);

        self.grid = space;
        return self;
    }

    fn deinit(self: *GridSamples, allocator: mem.Allocator) void {
        allocator.free(self.vals);
        allocator.destroy(self);
    }
};

test "test" {
    const allocator = testing.allocator;

    const interval: Interval = try Interval.init(0.0, 0.7);
    const chebyshev_grid: *ChebyshevGrid = try ChebyshevGrid.init(allocator, 3);
    defer chebyshev_grid.deinit(allocator);

    debug.print("{any}\n", .{chebyshev_grid.grid});

    const space_grid: *SpaceGrid = try SpaceGrid.initUninitialized(allocator, 3);
    defer space_grid.deinit(allocator);
    try chebyshev_grid.chebyshevToSpace(space_grid, interval);

    debug.print("{any}\n", .{space_grid.grid});

    const source_sampling: *GridSamples = try GridSamples.init(allocator, space_grid);
    defer source_sampling.deinit(allocator);
}

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
const math = std.math;
const debug = std.debug;
const testing = std.testing;
