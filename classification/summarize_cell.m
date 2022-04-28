function summarize_cell(charVector)

charVector = categorical(charVector);
num_samples = countcats(charVector);
[C,new_order] = sort(num_samples,'descend');
charVector = reordercats(charVector,new_order);
num_samples = countcats(charVector);
sample_type = unique(charVector);
summary(charVector);

return