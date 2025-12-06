function x_comb = comb(xf_m, xf_gm, th)
x_comb = xf_gm; pos_gm = find(x_comb > th);
for i = 1:length(pos_gm)
gm_pos = pos_gm(i); left_bound = gm_pos; right_bound = gm_pos;

while left_bound > 1 && xf_m(left_bound - 1) > th
    left_bound = left_bound - 1;
end
while right_bound < length(xf_m) && xf_m(right_bound + 1) > th
    right_bound = right_bound + 1;
end
x_comb(left_bound:right_bound) = xf_m(left_bound:right_bound);
i=right_bound;
end
end
