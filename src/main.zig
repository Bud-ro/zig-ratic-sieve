const std = @import("std");

pub fn main() !void {
    var gpa: std.heap.GeneralPurposeAllocator(.{}) = .init;
    const alloc = gpa.allocator();

    var arg_it = try std.process.argsWithAllocator(alloc);
    defer arg_it.deinit();

    // Assume first argument is the number to be factored
    const n = blk: {
        _ = arg_it.next().?;
        const num_string = arg_it.next().?;
        const n = try std.fmt.parseInt(u64, num_string, 10);
        break :blk n;
    };

    const stdout_file = std.io.getStdOut().writer();
    var bw = std.io.bufferedWriter(stdout_file);
    const stdout = bw.writer();

    try stdout.print("Factoring {}...\n", .{n});
    const factors = try factor(n, alloc);
    defer alloc.free(factors);

    try stdout.print("The factors are:", .{});
    for (factors) |q| {
        try stdout.print(" {}", .{q});
    }
    try stdout.print("\n", .{});

    try bw.flush();
}

/// Returns a slice of sorted factors that the caller owns
fn factor(n: u64, allocator: std.mem.Allocator) ![]u64 {
    // Fermat Factorization to begin
    var i: u64 = 0;
    while (i < n) : (i += 1) {
        const sqrt_candidate: u64 = std.math.sqrt(n + i * i);
        if (sqrt_candidate * sqrt_candidate == n + i * i) {
            const perfect_square = sqrt_candidate;

            var factors = try allocator.alloc(u64, 2);
            factors[0] = perfect_square - i;
            factors[1] = perfect_square + i;
            return factors;
        }
    }

    return error.CouldNotFactor;
}

test "Can factor small numbers" {
    const result = try factor(10, std.testing.allocator);
    try std.testing.expectEqualSlices(u64, &[_]u64{ 2, 5 }, result);
    std.testing.allocator.free(result);
}

test "Can factor smallish numbers" {
    const result = try factor(8051, std.testing.allocator);
    try std.testing.expectEqualSlices(u64, &[_]u64{ 83, 97 }, result);
    std.testing.allocator.free(result);
}

// This is noticeably slow, so we already need to ditch Fermat factorization
test "Can factor medium numbers" {
    const result = try factor(1299709 * 15485863, std.testing.allocator);
    try std.testing.expectEqualSlices(u64, &[_]u64{ 1299709, 15485863 }, result);
    std.testing.allocator.free(result);
}
