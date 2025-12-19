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

    fn chebyshevToSpace(self: *const ChebyshevGrid, space: *SpaceGrid, interval: Interval) !void {
        if (space.grid.len != self.grid.len) return GridError.LengthMismatch;
        for (space.grid, self.grid) |*dest, src| dest.* = interval.chebyshevToSpace(src);
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

    fn spaceToChebyshev(self: *const SpaceGrid, chebyshev: *ChebyshevGrid, interval: Interval) !void {
        if (chebyshev.grid.len != self.grid.len) return GridError.LengthMismatch;
        for (chebyshev.grid, self.grid) |*dest, src| dest.* = interval.spaceToChebyshev(src);
    }
};

/// We'll likely have lots of sampled things (u( t[l] ), G( t[l] ), s( t[l] ), etc.)
const GridSamples = struct {
    vals: []f64, // g( y[l] ) = f( t[l] ) values

    fn initUninitialized(allocator: mem.Allocator, order: usize) !*GridSamples {
        if (order < 2) return GridError.InvalidOrder;

        const self: *GridSamples = try allocator.create(GridSamples);
        errdefer allocator.destroy(self);

        self.vals = try allocator.alloc(f64, order);
        errdefer allocator.free(self.vals);

        return self;
    }

    fn deinit(self: *GridSamples, allocator: mem.Allocator) void {
        allocator.free(self.vals);
        allocator.destroy(self);
    }

    fn fillWith(self: *GridSamples, space: *const SpaceGrid, f: *const fn (t: f64) f64) !void {
        if (self.vals.len != space.grid.len) return GridError.LengthMismatch;
        for (self.vals, space.grid) |*v, t| v.* = f(t);
    }
};

const Backend = enum { chebyshev, newton };

const Coefficients = struct {
    // for Chebyshev: c[0..M-1]
    // for Newton:    a[0..M-1]
    data: []f64,

    // explicit: only used by Newton-stability transform (length-4 mapping)
    scale: f64,

    fn deinit(self: *Coefficients, allocator: mem.Allocator) void {
        allocator.free(self.data);
        allocator.destroy(self);
    }
};

fn Polynomial(comptime backend: Backend) type {
    return struct {
        fn build(allocator: mem.Allocator, chebyshev: *const ChebyshevGrid, samples: *const GridSamples, interval: Interval) !*Coefficients {
            _ = interval;

            if (chebyshev.grid.len != samples.vals.len) return GridError.LengthMismatch;

            const data: []f64 = try allocator.alloc(f64, samples.vals.len);
            errdefer allocator.free(data);

            switch (backend) {
                .chebyshev => @compileError("Not supported yet."),
                .newton => {
                    @memcpy(data, samples.vals);
                    const x: []f64 = chebyshev.grid;

                    for (1..data.len) |j| {
                        var i: usize = data.len - 1;
                        while (true) : (i -= 1) {
                            data[i] = (data[i] - data[i - 1]) / (x[i] - x[i - j]);
                            if (i == j) break;
                        }
                    }
                },
            }

            const coeff: *Coefficients = try allocator.create(Coefficients);
            errdefer allocator.destroy(coeff);

            coeff.data = data;

            // scale is reserved for Newton domain scaling; currently always 1.0
            coeff.scale = switch (backend) {
                .chebyshev, .newton => 1.0,
            };

            return coeff;
        }

        // evaluate at normalized y (can be slightly outside [-1,1] if caller wants extrapolation)
        fn eval(ynew: f64, coeffs: *const Coefficients, chebyshev: *const ChebyshevGrid) f64 {
            const y: []f64 = chebyshev.grid;
            const n: usize = y.len;

            switch (backend) {
                .chebyshev => @compileError("Not supported yet."),
                .newton => {
                    const a: []f64 = coeffs.data;
                    debug.assert(n == a.len); // user-unreachable

                    var i: usize = n - 1;
                    var p: f64 = a[i];

                    while (0 < i) {
                        i -= 1;
                        p = p * (ynew - y[i]) + a[i];
                    }

                    return p;
                },
            }
        }
    };
}

test "Newton polynomial interpolates samples at Chebyshev grids" {
    const allocator = testing.allocator;
    const order: comptime_int = 5;

    const interval: Interval = try Interval.init(0.0, 0.7);

    const chebyshev_grid: *ChebyshevGrid = try ChebyshevGrid.init(allocator, order);
    defer chebyshev_grid.deinit(allocator);

    // check Chebyshev endpoints
    try testing.expectApproxEqAbs(-1.0, chebyshev_grid.grid[0], 0.0);
    try testing.expectApproxEqAbs(1.0, chebyshev_grid.grid[order - 1], 0.0);

    const space_grid: *SpaceGrid = try SpaceGrid.initUninitialized(allocator, order);
    defer space_grid.deinit(allocator);
    try chebyshev_grid.chebyshevToSpace(space_grid, interval);

    // check mapped-space endpoints
    try testing.expectApproxEqAbs(interval.t0, space_grid.grid[0], 1e-15);
    try testing.expectApproxEqAbs(interval.t0 + interval.dt, space_grid.grid[order - 1], 1e-15);

    // check monotonicity
    for (space_grid.grid[0 .. order - 1], space_grid.grid[1..]) |t0, t1| {
        try testing.expect(t0 < t1);
    }

    const grid_samples: *GridSamples = try GridSamples.initUninitialized(allocator, order);
    defer grid_samples.deinit(allocator);
    try grid_samples.fillWith(space_grid, simpleSourceForTest);

    const NewtonPoly: type = Polynomial(.newton);
    const coeffs: *Coefficients = try NewtonPoly.build(allocator, chebyshev_grid, grid_samples, interval);
    defer coeffs.deinit(allocator);

    // Check p( y[l] ) == f( t[l] ) at all grid points
    for (chebyshev_grid.grid, grid_samples.vals) |y, f| {
        try testing.expectApproxEqAbs(f, NewtonPoly.eval(y, coeffs, chebyshev_grid), 1e-15);
    }
}

// Simple source: s(t) = 2 - 3t
fn simpleSourceForTest(t: f64) f64 {
    return 2.0 - 3.0 * t;
}

const std = @import("std");
const mem = std.mem;
const math = std.math;
const debug = std.debug;
const testing = std.testing;
