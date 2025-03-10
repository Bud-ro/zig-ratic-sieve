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

    // The maximum prime allowed as a factor.
    const smoothness_bound = 19;
    const factor_base = try generate_factor_base(alloc, smoothness_bound);
    defer alloc.free(factor_base);

    // Find smooth number candidates
    const smooth_numbers = try sieve(n, factor_base, 10000, alloc);
    defer alloc.free(smooth_numbers);

    const stdout_file = std.io.getStdOut().writer();
    var bw = std.io.bufferedWriter(stdout_file);
    const stdout = bw.writer();

    try stdout.print("Smoothness bound: {}\n", .{smoothness_bound});
    try stdout.print("Factor base: {d}\n", .{factor_base});
    try stdout.print("Smooth numbers: {d}\n", .{smooth_numbers});

    try stdout.print("Factoring {}...\n", .{n});
    try bw.flush();

    const factors = try factor(n, alloc);
    defer alloc.free(factors);

    try stdout.print("The factors are:", .{});
    for (factors) |q| {
        try stdout.print(" {}", .{q});
    }
    try stdout.print("\n", .{});

    try bw.flush();
}

/// Returns a caller-owned slice of primes for the factor base
fn generate_factor_base(gpa: std.mem.Allocator, smoothness_bound: u64) ![]u64 {
    var factor_base: std.ArrayList(u64) = .init(gpa);
    defer factor_base.deinit();

    // Find primes up to smoothness bound with brute force:
    var i: u64 = 2;
    while (i <= smoothness_bound) : (i += 1) {
        if (is_prime_brute_force(i)) {
            try factor_base.append(i);
        }
    }

    return factor_base.toOwnedSlice();
}

// TODO: Eliminate this when we get around to sieving (and large smoothness bounds)
fn is_prime_brute_force(n: u64) bool {
    var i: u64 = 2;
    while (i < std.math.sqrt(n) + 1) : (i += 1) {
        if (n % i == 0) {
            return false;
        }
    }

    return true;
}

/// Uses the polynomial Q(x) = x^2 - n to sieve for smooth numbers around sqrt(n)
fn sieve(n: u64, factor_base: []const u64, sieve_size: u64, gpa: std.mem.Allocator) ![]u64 {
    var sieve_array = try gpa.alloc(i32, sieve_size);
    defer gpa.free(sieve_array);
    @memset(sieve_array, 0);

    // Find quadratic residues for each prime in the factor base
    for (factor_base) |p| {
        // Find roots where x^2 ≡ n (mod p)
        var x: u64 = 0;
        while (x < p) : (x += 1) {
            // Prevent overflow by avoiding x*x % p
            if (((x % p) * (x % p)) % p == n % p) {
                break;
            }
        }

        if (x == p) continue; // No solution found, skip prime

        // Mark sieve positions corresponding to solutions
        var k: u64 = x;
        while (k < sieve_size) : (k += p) {
            sieve_array[k] += 1; // Increment count of factor base divisors
        }
    }

    // Collect indices with enough factors (potential smooth numbers)
    var smooth_numbers = std.ArrayList(u64).init(gpa);
    defer smooth_numbers.deinit();

    for (sieve_array, 0..) |count, index| {
        if (count > 3) { // Arbitrary threshold
            try smooth_numbers.append(index);
        }
    }

    return try smooth_numbers.toOwnedSlice();
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
