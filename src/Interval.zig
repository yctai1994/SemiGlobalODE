t0: f64,
dt: f64,

const Self: type = @This();

const ErrorSet = error{InvalidInterval};

pub fn init(t0: f64, dt: f64) ErrorSet!Self {
    if (dt <= 0.0) return ErrorSet.InvalidInterval;
    return .{ .t0 = t0, .dt = dt };
}

// map t -> y in [-1,1]
inline fn spaceToChebyshev(self: *const Self, t: f64) f64 {
    return (2.0 / self.dt) * (t - self.t0) - 1;
}

// map y -> t
pub inline fn chebyshevToSpace(self: *const Self, y: f64) f64 {
    return self.t0 + 0.5 * self.dt * (1 + y);
}
