SCALE = 10;
edgefactor = 16;
NBFS = 64;
edges_fname = "ijw.dat";
keys_fname = "keys.dat";
hist_fname = "hist.dat";
parent_fname = "parent.dat";
nedges_fname = "n_edges.dat";

rand ("seed", 103);

ijw = kronecker_generator (SCALE, edgefactor);

num_vertices = 2^SCALE;
num_edges = num_vertices * edgefactor;

f = fopen(edges_fname, "w");
fprintf(f, "%d %d\n", num_vertices, num_edges);
for i=1:1:num_edges
  fprintf(f, "%d %d %f\n", ijw(1,i)-1, ijw(2,i)-1, ijw(3,i));
end
fclose(f);

tic;
G = kernel_1 (ijw);
kernel_1_time = toc;

N = size (G, 1);
coldeg = full (spstats (G));
search_key = randperm (N);

search_key(coldeg(search_key) == 0) = [];
if length (search_key) > NBFS,
  search_key = search_key(1:NBFS);
else
  NBFS = length (search_key);
end
search_key = search_key - 1;

f = fopen(keys_fname, "w");
for i=1:1:length(search_key)
  fprintf(f, "%d\n", search_key(i));
end
fclose(f);

kernel_2_time = Inf * ones (NBFS, 1);
kernel_2_nedge = zeros (NBFS, 1);
kernel_3_time = Inf * ones (NBFS, 1);
kernel_3_nedge = zeros (NBFS, 1);

indeg = histc (ijw(:), 1:N); % For computing the number of edges

f = fopen(hist_fname, "w");
for i=1:1:N
  fprintf(f, "%d\n", indeg(i));
end
fclose(f);

f = fopen(parent_fname, "w");
f2 = fopen(nedges_fname, "w");

for k = 1:NBFS,
  tic;
  parent = kernel_2 (G, search_key(k), 2^SCALE);
  for p=1:1:length(parent)
    fprintf(f, "%d", parent(p))
    if (length(parent) == p) fprintf(f, "\n")
    else
      fprintf(f, " ");
    end
  end
  kernel_2_time(k) = toc;
  err = validate (parent, ijw, search_key(k), 0, false);
  if err <= 0,
    error (sprintf ("BFS %d from search key %d failed to validate: %d",
      k, search_key(k), err));
  end
  kernel_2_nedge(k) = sum (indeg(parent >= 0))/2; % Volume/2
  fprintf(f2, "%.1f %.4f\n", kernel_2_nedge(k), kernel_2_time(k))

  tic;
  [parent, d] = kernel_3 (G, search_key(k));
  kernel_3_time(k) = toc;
  err = validate (parent, ijw, search_key (k), d, true);
  if err <= 0,
    error (sprintf ("SSSP %d from search key %d failed to validate: %d",
      k, search_key(k), err));
  end
  kernel_3_nedge(k) = sum (indeg(parent >= 0))/2; % Volume/2
end

fclose(f);
fclose(f2);
output (SCALE, NBFS, NBFS, kernel_1_time, kernel_2_time, kernel_2_nedge, kernel_3_time, kernel_3_nedge);
