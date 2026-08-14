function d = correlation_matrix_distance(q1,Q1,q2,Q2)
% Correlation matrix distance from 
% https://www.researchgate.net/publication/4194743_Correlation_Matrix_Distance_a_Meaningful_Measure_for_Evaluation_of_Non-Stationary_MIMO_Channels
% d = 0 if matrices are the same, is 1 at most if matrices are different

if isequal(q1,q2)
    d = 1 - trace(Q1*Q2)/(norm(Q1,'fro')*norm(Q2,'fro'));
else
    d = 1;
end
end