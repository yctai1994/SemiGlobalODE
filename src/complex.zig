pub const ComplexF64 = struct {
    re: f64,
    im: f64,

    pub const zero: ComplexF64 = .{ .re = 0.0, .im = 0.0 };
    pub const one: ComplexF64 = .{ .re = 1.0, .im = 0.0 };

    pub inline fn show(z: *const ComplexF64) void {
        if (z.im < 0.0) {
            debug.print("{d: >9.5} - {d: >9.5}im", .{ z.re, @abs(z.im) });
        } else {
            debug.print("{d: >9.5} + {d: >9.5}im", .{ z.re, z.im });
        }
    }
};

pub inline fn complex(re: f64, im: f64) ComplexF64 {
    return .{ .re = re, .im = im };
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

pub fn add(z1: anytype, z2: anytype) ComplexF64 {
    const T1: type = comptime @TypeOf(z1);
    const T2: type = comptime @TypeOf(z2);

    switch (T1) {
        f64, comptime_float => {
            return switch (T2) {
                f64, comptime_float => complex(z1 + z2, 0.0),
                ComplexF64 => complex(z1 + z2.re, z2.im),
                else => unreachable,
            };
        },
        ComplexF64 => {
            return switch (T2) {
                f64, comptime_float => complex(z1.re + z2, z1.im),
                ComplexF64 => complex(z1.re + z2.re, z1.im + z2.im),
                else => unreachable,
            };
        },
        else => unreachable,
    }
}

test "ComplexF64 Addition" {
    const x: ComplexF64 = complex(3.0, 2.0);
    const y: ComplexF64 = complex(2.0, 1.0);

    try testing.expectEqual(complex(5.0, 0.0), add(3.0, 2.0));
    try testing.expectEqual(complex(5.0, 2.0), add(x, 2.0));
    try testing.expectEqual(complex(5.0, 1.0), add(3.0, y));
    try testing.expectEqual(complex(5.0, 3.0), add(x, y));
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

pub fn sub(z1: anytype, z2: anytype) ComplexF64 {
    const T1: type = comptime @TypeOf(z1);
    const T2: type = comptime @TypeOf(z2);

    switch (T1) {
        f64, comptime_float => {
            return switch (T2) {
                f64, comptime_float => complex(z1 - z2, 0.0),
                ComplexF64 => complex(z1 - z2.re, -z2.im),
                else => unreachable,
            };
        },
        ComplexF64 => {
            return switch (T2) {
                f64, comptime_float => complex(z1.re - z2, z1.im),
                ComplexF64 => complex(z1.re - z2.re, z1.im - z2.im),
                else => unreachable,
            };
        },
        else => unreachable,
    }
}

test "ComplexF64 Subtraction" {
    const x: ComplexF64 = complex(3.0, 2.0);
    const y: ComplexF64 = complex(2.0, -1.0);

    try testing.expectEqual(complex(1.0, 0.0), sub(3.0, 2.0));
    try testing.expectEqual(complex(1.0, 2.0), sub(x, 2.0));
    try testing.expectEqual(complex(1.0, 1.0), sub(3.0, y));
    try testing.expectEqual(complex(1.0, 3.0), sub(x, y));
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

pub fn mul(z1: anytype, z2: anytype) ComplexF64 {
    const T1: type = comptime @TypeOf(z1);
    const T2: type = comptime @TypeOf(z2);

    switch (T1) {
        f64, comptime_float => {
            return switch (T2) {
                f64, comptime_float => complex(z1 * z2, 0.0),
                ComplexF64 => complex(z1 * z2.re, z1 * z2.im),
                else => unreachable,
            };
        },
        ComplexF64 => {
            return switch (T2) {
                f64, comptime_float => complex(z1.re * z2, z1.im * z2),
                ComplexF64 => complex(
                    z1.re * z2.re - z1.im * z2.im,
                    z1.re * z2.im + z1.im * z2.re,
                ),
                else => unreachable,
            };
        },
        else => unreachable,
    }
}

test "ComplexF64 Multiplication" {
    const x: ComplexF64 = complex(3.0, 2.0);
    const y: ComplexF64 = complex(2.0, 1.0);

    try testing.expectEqual(complex(6.0, 0.0), mul(3.0, 2.0));
    try testing.expectEqual(complex(6.0, 4.0), mul(x, 2.0));
    try testing.expectEqual(complex(6.0, 3.0), mul(3.0, y));
    try testing.expectEqual(complex(4.0, 7.0), mul(x, y));
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

pub fn div(z1: anytype, z2: anytype) ComplexF64 {
    const T1: type = comptime @TypeOf(z1);
    const T2: type = comptime @TypeOf(z2);

    switch (T1) {
        f64, comptime_float => {
            return switch (T2) {
                f64, comptime_float => return complex(z1 / z2, 0.0),
                ComplexF64 => {
                    const de: f64 = sqr(z2.re) + sqr(z2.im);
                    return complex(z1 * z2.re / de, -z1 * z2.im / de);
                },
                else => unreachable,
            };
        },
        ComplexF64 => {
            switch (T2) {
                f64, comptime_float => return complex(z1.re / z2, z1.im / z2),
                ComplexF64 => {
                    const de: f64 = sqr(z2.re) + sqr(z2.im);
                    return complex(
                        (z1.re * z2.re + z1.im * z2.im) / de,
                        (z1.im * z2.re - z1.re * z2.im) / de,
                    );
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

test "ComplexF64 Division" {
    const x: ComplexF64 = complex(3.0, 2.0);
    const y: ComplexF64 = complex(2.0, -1.0);

    try testing.expectEqual(complex(1.5, 0.0), div(3.0, 2.0));
    try testing.expectEqual(complex(1.5, 1.0), div(x, 2.0));
    try testing.expectEqual(complex(1.2, 0.6), div(3.0, y));
    try testing.expectEqual(complex(0.8, 1.4), div(x, y));
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

pub inline fn abs(z: ComplexF64) f64 {
    return @sqrt(abs2(z));
}

pub inline fn cis(x: f64) ComplexF64 {
    return complex(@cos(x), @sin(x));
}

pub inline fn conj(z: ComplexF64) ComplexF64 {
    return complex(z.re, -z.im);
}

pub fn abs2(z: ComplexF64) f64 {
    const x: f64 = z.re;
    if (math.isNan(x)) return x;
    const y: f64 = z.im;
    if (math.isNan(y)) return y;
    // general case
    const X: f64 = @abs(x);
    const Y: f64 = @abs(y);
    const u: f64 = @max(X, Y);
    const v: f64 = @min(X, Y);
    if (v == 0.0) return sqr(u);
    return sqr(u) * (1.0 + sqr(v / u));
}

pub inline fn angle(z: ComplexF64) f64 {
    return math.atan2(z.im, z.re);
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

inline fn sqr(x: f64) f64 {
    return x * x;
}

const std = @import("std");
const math = std.math;
const debug = std.debug;
const testing = std.testing;
