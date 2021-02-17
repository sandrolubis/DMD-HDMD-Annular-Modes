function W = my_compress(V, q)

    m = size(V, 2);
    n = size(V, 1)/q;

    W = zeros(n, m);

%% My way    
    
% for(j = 1:m)
%     for(i = 1:n)
%         entry = 0;
%         for(k = 1:q)
%             entry = entry + V((i-1)*q+k, j);
%         end
%         W(i, j) = entry/q;
%     end
% end

%% Mezic way

% for(j = 1:m)
%     for(i = 1:n)
%         entry = 0;
%         for(k = 1:q)
%             entry = entry + V((k-1)*n+i, j);
%         end
%         W(i, j) = entry/q;
%     end
% end


for(j = 1:m)
    W(:, j) = V(end-n+1:end, j);
end


return