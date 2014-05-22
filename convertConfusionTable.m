% convert a confusion matrix into a confusion table of true/false positives/negatives
function conf_table = convertConfusionTable(cfm, label)
% Returns a confusion table in the following format:
% [true positives, false positives, true+false positives;
% false negatives, true negatives, false+true negatives;
% true+false positives, false+true positives, true positives+false positives+false negatives+true negatives]
% for the given label index in the confusion matrix.

predicted = cfm(:, label);
actual = cfm(label, :);
% for ii = 1:length(cfm)
%     actual(1, ii) = cfm(ii, label);
% end
true_pos  = predicted(label);
false_pos = sum(actual) - true_pos;
false_neg = sum(predicted) - true_pos;
total     = sum(sum(cfm, 2)); %sum column values
true_neg  = total - true_pos - false_pos - false_neg;

conf_table = [true_pos, false_pos; false_neg, true_neg];
end