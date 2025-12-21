/// "Space" means the physical independent variable (t, x, etc.)
grid: []f64,

const Self: type = @This();

const ErrorSet = error{ InvalidOrder, LengthMismatch } || mem.Allocator.Error;

pub fn initUninitialized(allocator: mem.Allocator, order: usize) ErrorSet!Self {
    if (order < 2) return ErrorSet.InvalidOrder;
    return .{ .grid = try allocator.alloc(f64, order) };
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

fn spaceToChebyshev(self: *const Self, chebyshev: *const ChebyshevGrid, interval: Interval) ErrorSet!void {
    if (chebyshev.grid.len != self.grid.len) return ErrorSet.LengthMismatch;
    for (chebyshev.grid, self.grid) |*dest, src| dest.* = interval.spaceToChebyshev(src);
}

const std = @import("std");
const mem = std.mem;
const Interval = @import("./Interval.zig");
const ChebyshevGrid = @import("./ChebyshevGrid.zig");
