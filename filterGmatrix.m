function [G, G_ind, related, n_genes_KO] = filterGmatrix(G, G_ind, maxlength, G_time, G_rxns, model_name)

Gind2genes_genes = unique([G_ind{:}]);
Gind2genes_mat = spalloc(length(G_ind),length(Gind2genes_genes), sum(cellfun(@length,G_ind)));
for i = 1:length(G_ind)
    Gind2genes_mat(i,ismember(Gind2genes_genes, G_ind{i})) = 1;
end
tabulate(sum(Gind2genes_mat,2)<=maxlength)
pos = sum(Gind2genes_mat,2)<=maxlength;
G = G(pos,:);
G_ind = G_ind(pos,1);

n_genes_KO = cellfun(@length, G_ind);
[n_genes_KO, ind] = sort(n_genes_KO, 'ascend');
G_ind = G_ind(ind);
G = G(ind, :);
n_G_ind = length(G_ind);

related = zeros(0,2);

% generate matrix that relates genes and G_ind
Gind2genes_genes = unique([G_ind{:}]);
Gind2genes_mat = spalloc(length(G_ind),length(Gind2genes_genes), sum(cellfun(@length,G_ind)));
for i = 1:length(G_ind)
    Gind2genes_mat(i,ismember(Gind2genes_genes, G_ind{i})) = 1;
end
% use matrix to search G_inds that contains lower order G_inds
for i = 1:n_G_ind
    act_G_ind = G_ind{i};
    n_act_G_ind = length(act_G_ind);
    pos = find(n_genes_KO > n_act_G_ind);
    % if a G_ind contains a smaller order one, all columns for these genes
    % should be one
    pos = pos(mean(Gind2genes_mat(pos,ismember(Gind2genes_genes,act_G_ind)),2)==1);
    
    for j = 1:length(pos)
        % Increase the G matrix
        G(pos(j), :) = G(pos(j), :) + G(i, :);
        
        % Add the relationships between G_inds
        related(end+1,:) = [pos(j), i];
    end
end
G = double(G>0);

if size(related)>0
    un_related = unique(related(:, 1));
    n_un_related = length(un_related);
    for i = 1:n_un_related
        act_KO = un_related(i);
        ind = related(:, 1) == act_KO;
        act_related = related(ind, 2);
        all_genes = [G_ind{act_related}];
        un_all_genes = unique(all_genes);
        n_un_all_genes = length(un_all_genes);
        n_genes_KO(act_KO) = n_genes_KO(act_KO)-n_un_all_genes;
    end
else
    related = NaN;
end

final_filename = [pwd filesep 'G_' model_name '.mat'];

save(final_filename, 'G', 'G_ind', 'G_rxns', 'related', 'n_genes_KO', 'G_time');

end
