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

    const exponent_matrix = try build_exponent_matrix(smooth_numbers, factor_base, alloc);
    defer {
        for (exponent_matrix) |row| {
            alloc.free(row);
        }
        alloc.free(exponent_matrix);
    }

    gaussian_elimination_mod2(exponent_matrix);
    const dependent_rows = try find_dependent_rows(exponent_matrix, alloc);

    const stdout_file = std.io.getStdOut().writer();
    var bw = std.io.bufferedWriter(stdout_file);
    const stdout = bw.writer();

    try stdout.print("Smoothness bound: {}\n", .{smoothness_bound});
    try stdout.print("Factor base: {d}\n", .{factor_base});
    try stdout.print("Smooth numbers: {d}\n", .{smooth_numbers});
    try stdout.print("Dependent rows: {d}\n", .{dependent_rows});

    try stdout.print("Factoring {}...\n", .{n});
    try bw.flush();

    // Attempt to factor using the matrix
    if (dependent_rows.len > 0) {
        const factors = try extract_factors(n, smooth_numbers, factor_base, dependent_rows, alloc);
        if (factors.len > 0) {
            std.debug.print("Nontrivial factors found: {d}\n", .{factors});
        } else {
            std.debug.print("No factor found using this dependency set.\n", .{});
        }
    } else {
        std.debug.print("No dependent rows found, try increasing smoothness bound.\n", .{});
    }

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
        if (count > 2) { // Arbitrary threshold
            try smooth_numbers.append(index);
        }
    }

    return try smooth_numbers.toOwnedSlice();
}

/// Returns a binary matrix of size factor_base.len x factor_base.len
fn build_exponent_matrix(smooth_numbers: []const u64, factor_base: []const u64, allocator: std.mem.Allocator) ![][]u8 {
    var matrix = try allocator.alloc([]u8, smooth_numbers.len);

    for (smooth_numbers, 0..) |num, i| {
        var row = try allocator.alloc(u8, factor_base.len);
        @memset(row, 0);

        var remainder = num;
        for (factor_base, 0..) |p, j| {
            var exponent: u8 = 0;
            while (remainder % p == 0) {
                remainder /= p;
                exponent += 1;
            }
            row[j] = exponent % 2; // Store only parity of exponents
        }

        matrix[i] = row;
    }

    return matrix;
}

/// Performs Gaussian elimination on our exponent mod 2 matrix
fn gaussian_elimination_mod2(matrix: [][]u8) void {
    const row_count = matrix.len;
    const col_count = matrix[0].len;

    var row: usize = 0;
    for (0..col_count) |col| {
        var pivot_row: ?usize = null;

        // Find a row with a leading 1 in this column
        for (row..row_count) |r| {
            if (matrix[r][col] == 1) {
                pivot_row = r;
                break;
            }
        }

        if (pivot_row == null) continue; // No pivot, skip column

        // Swap pivot row into position
        std.mem.swap([]u8, &matrix[row], &matrix[pivot_row.?]);

        // Eliminate 1s in this column below the pivot
        for (row + 1..row_count) |r| {
            if (matrix[r][col] == 1) {
                for (0..col_count) |c| {
                    matrix[r][c] ^= matrix[row][c]; // XOR row with pivot row
                }
            }
        }

        row += 1; // Move to next row
    }
}

fn find_dependent_rows(matrix: [][]u8, gpa: std.mem.Allocator) ![]usize {
    var dependent_rows = std.ArrayList(usize).init(gpa);
    defer dependent_rows.deinit();

    outer: for (matrix, 0..) |row, i| {
        // If all zeros, it's dependent
        for (row) |elem| {
            if (elem != 0) {
                continue :outer;
            }
        }
        try dependent_rows.append(i);
    }

    return dependent_rows.toOwnedSlice();
}

/// Extracts nontrivial factors from the exponent matrix dependencies
fn extract_factors(n: u64, smooth_numbers: []const u64, factor_base: []const u64, dependent_rows: []const usize, allocator: std.mem.Allocator) ![]u64 {
    var x: u64 = 1;
    var y: u64 = 1;
    var factors = std.ArrayList(u64).init(allocator);
    defer factors.deinit();

    for (dependent_rows) |i| {
        var num = smooth_numbers[i];
        for (factor_base) |p| {
            var exponent: u8 = 0;
            while (num % p == 0) {
                num /= p;
                exponent += 1;
            }
            if (exponent % 2 == 1) {
                x *= p;
            }
        }
        y *= smooth_numbers[i];
    }

    x = std.math.sqrt(x);
    y %= n;

    const factor1 = std.math.gcd(x + y, n);
    const factor2 = std.math.gcd(x - y, n);

    if (factor1 != 1 and factor1 != n) try factors.append(factor1);
    if (factor2 != 1 and factor2 != n) try factors.append(factor2);

    if (factors.items.len == 0) return error.CouldNotFactor;
    return factors.toOwnedSlice();
}

// test "Can factor small numbers" {
//     const result = try factor(10, std.testing.allocator);
//     try std.testing.expectEqualSlices(u64, &[_]u64{ 2, 5 }, result);
//     std.testing.allocator.free(result);
// }

// test "Can factor smallish numbers" {
//     const result = try factor(8051, std.testing.allocator);
//     try std.testing.expectEqualSlices(u64, &[_]u64{ 83, 97 }, result);
//     std.testing.allocator.free(result);
// }

// // This is noticeably slow, so we already need to ditch Fermat factorization
// test "Can factor medium numbers" {
//     const result = try factor(1299709 * 15485863, std.testing.allocator);
//     try std.testing.expectEqualSlices(u64, &[_]u64{ 1299709, 15485863 }, result);
//     std.testing.allocator.free(result);
// }
