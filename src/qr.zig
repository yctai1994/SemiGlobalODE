test "test" {
    const page = testing.allocator;

    const ArrComplexF64: ArrayComplexF64 = .{ .allocator = page };

    const N: comptime_int = 3;

    const R: [][]ComplexF64 = try ArrComplexF64.matrix(N, N);
    defer ArrComplexF64.free(R);

    inline for (R[0], .{ 1.0, -5.0, 4.0 }, .{ 2.0, 1.0, 0.0 }) |*p, re, im| {
        p.* = complex(re, im);
    }
    inline for (R[1], .{ 6.0, 16.0, -6.0 }, .{ 0.0, 7.0, 8.0 }) |*p, re, im| {
        p.* = complex(re, im);
    }
    inline for (R[2], .{ -4.0, 2.0, -4.0 }, .{ 0.0, 4.0, 1.0 }) |*p, re, im| {
        p.* = complex(re, im);
    }

    const Q: [][]ComplexF64 = try ArrComplexF64.matrix(N, N);
    defer ArrComplexF64.free(Q);

    for (Q, 0..) |row, i| {
        @memset(row, ComplexF64.zero);
        row[i] = ComplexF64.one;
    }

    debug.print("R = \n", .{});
    for (R) |row| {
        debug.print("[ ", .{});
        for (row) |val| {
            val.show();
            debug.print(", ", .{});
        }
        debug.print(" ]\n", .{});
    }

    debug.print("Q = \n", .{});
    for (Q) |row| {
        debug.print("[ ", .{});
        for (row) |val| {
            val.show();
            debug.print(", ", .{});
        }
        debug.print(" ]\n", .{});
    }

    for (0..N - 1) |k| {
        var scale: f64 = 0.0;
        for (k..N) |i| scale = @max(scale, Complex.abs(R[i][k]));
        for (k..N) |i| R[i][k] = Complex.div(R[i][k], scale);

        var sigma: f64 = 0.0;
        for (k..N) |i| sigma += Complex.abs2(R[i][k]);
        sigma = @sqrt(sigma);

        const phase: ComplexF64 = Complex.cis(Complex.angle(R[k][k]));
        const beta: f64 = 1.0 / (sigma * (sigma + Complex.abs(R[k][k])));

        R[k][k] = Complex.add(R[k][k], Complex.mul(sigma, phase));

        for (k + 1..N) |j| {
            var tmp: ComplexF64 = .zero;
            for (k..N) |i| {
                tmp = Complex.add(tmp, Complex.mul(Complex.conj(R[i][k]), R[i][j]));
            }
            tmp = Complex.mul(beta, tmp);

            for (k..N) |i| R[i][j] = Complex.sub(R[i][j], Complex.mul(tmp, R[i][k]));
        }

        for (0..N) |j| {
            var tmp: ComplexF64 = .zero;
            for (k..N) |i| {
                tmp = Complex.add(tmp, Complex.mul(Complex.conj(R[i][k]), Q[i][j]));
            }
            tmp = Complex.mul(beta, tmp);

            for (k..N) |i| Q[i][j] = Complex.sub(Q[i][j], Complex.mul(tmp, R[i][k]));
        }

        R[k][k] = Complex.mul(-scale * sigma, phase);
        for (k + 1..N) |i| R[i][k] = .zero;
    }

    debug.print("= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =\n", .{});

    debug.print("R = \n", .{});
    for (R) |row| {
        debug.print("[ ", .{});
        for (row) |val| {
            val.show();
            debug.print(", ", .{});
        }
        debug.print(" ]\n", .{});
    }

    debug.print("Q = \n", .{});
    for (Q) |row| {
        debug.print("[ ", .{});
        for (row) |val| {
            val.show();
            debug.print(", ", .{});
        }
        debug.print(" ]\n", .{});
    }
}

const std = @import("std");
const debug = std.debug;
const testing = std.testing;

const Complex = @import("./complex.zig");
const complex = Complex.complex;
const ComplexF64 = Complex.ComplexF64;
const ArrayComplexF64 = @import("./array.zig").Array(ComplexF64);
