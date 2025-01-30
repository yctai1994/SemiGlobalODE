fn qr(A: [][]ComplexF64, R: [][]ComplexF64, Qt: [][]ComplexF64, N: usize) void {
    for (0..N) |i| {
        @memcpy(R[i][0..N], A[i][0..N]);
        @memset(Qt[i][0..N], ComplexF64.zero);
        Qt[i][i] = ComplexF64.one;
    }

    for (0..N) |k| {
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
            for (k..N) |i| tmp = Complex.add(tmp, Complex.mul(Complex.conj(R[i][k]), R[i][j]));
            tmp = Complex.mul(beta, tmp);

            for (k..N) |i| R[i][j] = Complex.sub(R[i][j], Complex.mul(tmp, R[i][k]));
        }

        for (0..N) |j| {
            var tmp: ComplexF64 = .zero;
            for (k..N) |i| tmp = Complex.add(tmp, Complex.mul(Complex.conj(R[i][k]), Qt[i][j]));
            tmp = Complex.mul(beta, tmp);

            for (k..N) |i| Qt[i][j] = Complex.sub(Qt[i][j], Complex.mul(tmp, R[i][k]));
        }

        R[k][k] = Complex.mul(-scale * sigma, phase);
        for (k + 1..N) |i| R[i][k] = .zero;
    }
}

fn rq(A: [][]ComplexF64, R: [][]ComplexF64, Qt: [][]ComplexF64, N: usize) void {
    for (0..N) |i| {
        for (0..N) |j| {
            var tmp: ComplexF64 = .zero;
            for (i..N) |k| tmp = Complex.add(tmp, Complex.mul(R[i][k], Complex.conj(Qt[j][k])));
            A[i][j] = tmp;
        }
    }
}

fn eigen2(a: ComplexF64, b: ComplexF64, c: ComplexF64, d: ComplexF64, u: *ComplexF64, v: *ComplexF64) void {
    // ┌     ┐
    // │ a b │
    // │ c d │
    // └     ┘
    const m: ComplexF64 = Complex.mul(0.5, Complex.add(a, d));
    const p: ComplexF64 = Complex.sub(Complex.mul(a, d), Complex.mul(b, c));
    const y: ComplexF64 = Complex.sub(Complex.mul(m, m), p);

    const r: f64 = @sqrt(Complex.abs(y));
    const t: f64 = 0.5 * Complex.angle(y);

    const z: ComplexF64 = Complex.mul(r, Complex.cis(t));

    u.* = Complex.sub(m, z);
    v.* = Complex.add(m, z);
}

const Errors = error{HessenbergEigenFail};

pub fn eigen(w: []ComplexF64, A: [][]ComplexF64, R: [][]ComplexF64, Qt: [][]ComplexF64) !void {
    var n: usize = w.len;
    var p: usize = n - 1;
    var q: usize = n - 2;
    var c: usize = 0; // Iteration protection
    var s: ComplexF64 = undefined; // Wilkinson shift
    var u: ComplexF64 = undefined; // eigenvalue for Wilkinson shift
    var v: ComplexF64 = undefined; // eigenvalue for Wilkinson shift

    while (2 < n) : ({
        n -= 1;
        p -= 1;
        q -= 1;
        c = 0;
    }) {
        while (c < 30 and 0.0 < Complex.abs2(A[p][q])) : (c += 1) {
            eigen2(A[q][q], A[q][p], A[p][q], A[p][p], &u, &v);
            s = if (Complex.abs2(Complex.sub(u, A[p][p])) < Complex.abs2(Complex.sub(v, A[p][p]))) u else v;

            for (0..n) |i| A[i][i] = Complex.sub(A[i][i], s);

            qr(A, R, Qt, n);
            rq(A, R, Qt, n);

            for (0..n) |i| A[i][i] = Complex.add(A[i][i], s);
        } else {
            if (c == 30) return error.HessenbergEigenFail;
        }

        w[p] = A[p][p];
    }

    eigen2(A[q][q], A[q][p], A[p][q], A[p][p], &w[q], &w[p]);
}

test "test" {
    const page = testing.allocator;

    const ArrComplexF64: ArrayComplexF64 = .{ .allocator = page };

    const N: comptime_int = 5;

    const A: [][]ComplexF64 = try ArrComplexF64.matrix(N, N);
    defer ArrComplexF64.free(A);

    inline for (A[0], .{ 1.0, 5.0, 2.0, -3.0, 1.0 }, .{ 2.0, 1.0, 1.0, -1.0, 2.0 }) |*p, re, im| p.* = complex(re, im);
    inline for (A[1], .{ 3.0, 2.0, 6.0, -3.0, 1.0 }, .{ 1.0, 3.0, 7.0, -4.0, 5.0 }) |*p, re, im| p.* = complex(re, im);
    inline for (A[2], .{ 0.0, 1.0, 1.0, -1.0, 1.0 }, .{ 0.0, 5.0, 2.0, -3.0, 7.0 }) |*p, re, im| p.* = complex(re, im);
    inline for (A[3], .{ 0.0, 0.0, 7.0, -4.0, 3.0 }, .{ 0.0, 0.0, 6.0, -3.0, 1.0 }) |*p, re, im| p.* = complex(re, im);
    inline for (A[4], .{ 0.0, 0.0, 0.0, -2.0, 5.0 }, .{ 0.0, 0.0, 0.0, -2.0, 1.0 }) |*p, re, im| p.* = complex(re, im);

    const R: [][]ComplexF64 = try ArrComplexF64.matrix(N, N);
    defer ArrComplexF64.free(R);

    const Qt: [][]ComplexF64 = try ArrComplexF64.matrix(N, N);
    defer ArrComplexF64.free(Qt);

    const sol: []ComplexF64 = try ArrComplexF64.vector(N);
    defer ArrComplexF64.free(sol);

    try eigen(sol, A, R, Qt);

    const ans: [5]ComplexF64 = .{
        complex(-3.09127567631315, -4.956778631380267),
        complex(-2.1742167794942686, 0.5315174443854475),
        complex(-0.5674636931731067, 5.030327072493509),
        complex(5.841044038072055, 4.4142352393281525),
        complex(4.991912110908483, -0.019301124826846896),
    };

    for (ans, sol) |ans_val, sol_val| {
        try testing.expectApproxEqAbs(ans_val.re, sol_val.re, 1e-13);
        try testing.expectApproxEqAbs(ans_val.im, sol_val.im, 1e-13);
    }
}

const std = @import("std");
const debug = std.debug;
const testing = std.testing;

const Complex = @import("./complex.zig");
const complex = Complex.complex;
const ComplexF64 = Complex.ComplexF64;
const ArrayComplexF64 = @import("./array.zig").Array(ComplexF64);
