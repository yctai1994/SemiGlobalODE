const Backend = enum { chebyshev, newton };

const Coefficients = struct {
    // for Chebyshev: c[0..M-1]
    // for Newton:    a[0..M-1]
    data: []f64,

    // explicit: only used by Newton-stability transform (length-4 mapping)
    scale: f64,

    inline fn deinit(self: *const Coefficients, allocator: mem.Allocator) void {
        allocator.free(self.data);
    }
};

fn Polynomial(comptime backend: Backend) type {
    return struct {
        const Self: type = @This();

        const ErrorSet = error{ InvalidOrder, LengthMismatch } || mem.Allocator.Error;

        fn build(allocator: mem.Allocator, chebyshev: *const ChebyshevGrid, samples: *const GridSamples) ErrorSet!Coefficients {
            if (chebyshev.grid.len != samples.vals.len) return ErrorSet.LengthMismatch;

            const data: []f64 = try allocator.alloc(f64, samples.vals.len);
            errdefer allocator.free(data);

            switch (backend) {
                .chebyshev => @compileError("Not supported yet."),
                .newton => {
                    @memcpy(data, samples.vals);
                    const x: []const f64 = chebyshev.grid;

                    for (1..data.len) |j| {
                        var i: usize = data.len - 1;
                        while (true) : (i -= 1) {
                            data[i] = (data[i] - data[i - 1]) / (x[i] - x[i - j]);
                            if (i == j) break;
                        }
                    }
                },
            }

            return .{
                .data = data,
                .scale = switch (backend) {
                    .chebyshev, .newton => 1.0,
                },
            };
        }

        // evaluate at normalized y (can be slightly outside [-1,1] if caller wants extrapolation)
        fn eval(ynew: f64, coeffs: *const Coefficients, chebyshev: *const ChebyshevGrid) f64 {
            const y: []const f64 = chebyshev.grid;
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

test "Test on divided difference" {
    const allocator = testing.allocator;
    const order: comptime_int = 4;
    const slice: [order]f64 = .{ -2.0, -1.0, 3.0, 5.0 };

    const chebyshev_grid: ChebyshevGrid = try ChebyshevGrid.initFromValues(allocator, &slice);
    defer chebyshev_grid.deinit(allocator);

    const space_grid: SpaceGrid = try SpaceGrid.initFromValues(allocator, &slice);
    defer space_grid.deinit(allocator);

    const grid_samples: GridSamples = try GridSamples.initUninitialized(allocator, order);
    defer grid_samples.deinit(allocator);
    try grid_samples.fillWith(&space_grid, simplePolynomialForTest);

    const NewtonPoly: type = Polynomial(.newton);
    const coeffs: Coefficients = try NewtonPoly.build(allocator, &chebyshev_grid, &grid_samples);
    defer coeffs.deinit(allocator);

    inline for (.{ -3.0, -1.0, 2.0, 5.0 }, coeffs.data) |answer, result| {
        try testing.expectEqual(answer, result);
    }
}

test "Newton polynomial interpolates samples at Chebyshev grids" {
    const allocator = testing.allocator;
    const order: comptime_int = 5;

    const interval: Interval = try Interval.init(0.0, 0.7);

    const chebyshev_grid: ChebyshevGrid = try ChebyshevGrid.init(allocator, order);
    defer chebyshev_grid.deinit(allocator);

    // check Chebyshev endpoints
    try testing.expectApproxEqAbs(-1.0, chebyshev_grid.grid[0], 0.0);
    try testing.expectApproxEqAbs(1.0, chebyshev_grid.grid[order - 1], 0.0);

    const space_grid: SpaceGrid = try SpaceGrid.initUninitialized(allocator, order);
    defer space_grid.deinit(allocator);
    try chebyshev_grid.chebyshevToSpace(&space_grid, interval);

    // check mapped-space endpoints
    try testing.expectApproxEqAbs(interval.t0, space_grid.grid[0], 1e-15);
    try testing.expectApproxEqAbs(interval.t0 + interval.dt, space_grid.grid[order - 1], 1e-15);

    // check monotonicity
    for (space_grid.grid[0 .. order - 1], space_grid.grid[1..]) |t0, t1| {
        try testing.expect(t0 < t1);
    }

    const grid_samples: GridSamples = try GridSamples.initUninitialized(allocator, order);
    defer grid_samples.deinit(allocator);
    try grid_samples.fillWith(&space_grid, simplePolynomialForTest);

    const NewtonPoly: type = Polynomial(.newton);
    const coeffs: Coefficients = try NewtonPoly.build(allocator, &chebyshev_grid, &grid_samples);
    defer coeffs.deinit(allocator);

    // Check p( y[l] ) == f( t[l] ) at all grid points
    for (chebyshev_grid.grid, grid_samples.vals) |y, f| {
        try testing.expectApproxEqAbs(f, NewtonPoly.eval(y, &coeffs, &chebyshev_grid), 1e-15);
    }
}

fn simplePolynomialForTest(x: f64) f64 {
    const avec: [4]f64 = .{ -3.0, -1.0, 2.0, 5.0 };
    const xvec: [4]f64 = .{ -2.0, -1.0, 3.0, 5.0 };
    return avec[0] +
        avec[1] * (x - xvec[0]) +
        avec[2] * (x - xvec[0]) * (x - xvec[1]) +
        avec[3] * (x - xvec[0]) * (x - xvec[1]) * (x - xvec[2]);
}

const std = @import("std");
const mem = std.mem;
const debug = std.debug;
const testing = std.testing;

const Interval = @import("./Interval.zig");
const ChebyshevGrid = @import("./ChebyshevGrid.zig");
const SpaceGrid = @import("./SpaceGrid.zig");
const GridSamples = @import("./GridSamples.zig");
