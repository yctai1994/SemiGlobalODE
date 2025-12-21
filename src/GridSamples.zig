/// We'll likely have lots of sampled things (u( t[l] ), G( t[l] ), s( t[l] ), etc.)
vals: []f64, // g( y[l] ) = f( t[l] ) values

const Self: type = @This();

const ErrorSet = error{ InvalidOrder, LengthMismatch } || mem.Allocator.Error;

pub fn initUninitialized(allocator: mem.Allocator, order: usize) ErrorSet!Self {
    if (order < 2) return ErrorSet.InvalidOrder;
    return .{ .vals = try allocator.alloc(f64, order) };
}

pub inline fn deinit(self: *const Self, allocator: mem.Allocator) void {
    allocator.free(self.vals);
}

pub fn fillWith(self: *const Self, space: *const SpaceGrid, f: *const fn (t: f64) f64) ErrorSet!void {
    if (self.vals.len != space.grid.len) return ErrorSet.LengthMismatch;
    for (self.vals, space.grid) |*v, t| v.* = f(t);
}

const std = @import("std");
const mem = std.mem;
const SpaceGrid = @import("./SpaceGrid.zig");
