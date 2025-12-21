// "Chebyshev" means the dimensionless reference coordinate y ∈ [-1,1]
grid: []f64,

const Self: type = @This();

const ErrorSet = error{ InvalidOrder, LengthMismatch } || mem.Allocator.Error;

// boundary-including nodes y(l) = -cos( l*π / (M-1) )  (Eq. 76, J. Comp. Phys. 343, 368, 2017)
pub fn init(allocator: mem.Allocator, order: usize) ErrorSet!Self {
    if (order < 2) return ErrorSet.InvalidOrder;

    const grid: []f64 = try allocator.alloc(f64, order);
    errdefer allocator.free(grid);

    const de: f64 = @floatFromInt(order - 1);
    for (grid, 0..) |*ptr, ind| {
        const nu: f64 = @floatFromInt(ind);
        ptr.* = -@cos(math.pi * (nu / de));
    }

    return .{ .grid = grid };
}

pub fn initFromValues(allocator: mem.Allocator, values: []const f64) ErrorSet!Self {
    const order: usize = values.len;
    if (order < 2) return ErrorSet.InvalidOrder;

    const grid: []f64 = try allocator.alloc(f64, order);
    errdefer allocator.free(grid);

    @memcpy(grid, values);

    return .{ .grid = grid };
}

pub inline fn deinit(self: *const Self, allocator: mem.Allocator) void {
    allocator.free(self.grid);
}

pub fn chebyshevToSpace(self: *const Self, space: *const SpaceGrid, interval: Interval) ErrorSet!void {
    if (space.grid.len != self.grid.len) return ErrorSet.LengthMismatch;
    for (space.grid, self.grid) |*dest, src| dest.* = interval.chebyshevToSpace(src);
}

const std = @import("std");
const mem = std.mem;
const math = std.math;
const Interval = @import("./Interval.zig");
const SpaceGrid = @import("./SpaceGrid.zig");
