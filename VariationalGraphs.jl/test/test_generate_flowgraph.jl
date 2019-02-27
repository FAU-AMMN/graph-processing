f, c = generate_flowgraph(d, f0, W, F, g, alpha, "aniso")
part1, part2, flow = cutpursuit(d, f0, g, W, "aniso", alpha, 5);
