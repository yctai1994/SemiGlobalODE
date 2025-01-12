pub fn Array(comptime T: type) type {
    // already comptime scope
    const slice_al: comptime_int = @alignOf([]T);
    const child_al: comptime_int = @alignOf(T);
    const slice_sz: comptime_int = @sizeOf(usize) * 2;
    const child_sz: comptime_int = @sizeOf(T);

    return struct {
        allocator: std.mem.Allocator,

        const Self: type = @This();

        pub fn matrix(self: *const Self, nrow: usize, ncol: usize) std.mem.Allocator.Error![][]T {
            const buff: []u8 = try self.allocator.alloc(u8, nrow * ncol * child_sz + nrow * slice_sz);

            const mat: [][]T = blk: {
                const ptr: [*]align(slice_al) []T = @ptrCast(@alignCast(buff.ptr));
                break :blk ptr[0..nrow];
            };

            const chunk_sz: usize = ncol * child_sz;
            var padding: usize = nrow * slice_sz;

            for (mat) |*row| {
                row.* = blk: {
                    const ptr: [*]align(child_al) T = @ptrCast(@alignCast(buff.ptr + padding));
                    break :blk ptr[0..ncol];
                };
                padding += chunk_sz;
            }

            return mat;
        }

        pub fn vector(self: *const Self, n: usize) std.mem.Allocator.Error![]T {
            return try self.allocator.alloc(T, n);
        }

        pub fn free(self: *const Self, slice: anytype) void {
            const S: type = comptime @TypeOf(slice);

            switch (S) {
                [][]T => {
                    const ptr: [*]u8 = @ptrCast(@alignCast(slice.ptr));
                    const len: usize = blk: {
                        const nrow: usize = slice.len;
                        const ncol: usize = slice[0].len;
                        break :blk nrow * ncol * child_sz + nrow * slice_sz;
                    };

                    self.allocator.free(ptr[0..len]);
                },
                []T => {
                    self.allocator.free(slice);
                },
                else => @compileError("Invalid type: " ++ @typeName(T)),
            }

            return;
        }
    };
}

const std = @import("std");
