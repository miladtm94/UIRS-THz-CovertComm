function print_results(ATR,APC,AEE)
[mrows, ncols] = size([ATR;APC;AEE]);
outputstr = ['%1.' num2str(mrows) 'e  '];
template = ['%1.' num2str(mrows) 'e  '] ;
outputstr = template;
for i = 2:ncols
    outputstr = [outputstr template];
end
var_names = ["ATR"; "AFP"; "AEE"];
outputstr = [outputstr '\n'] ;
fprintf('\n');
fprintf(outputstr, [ATR;APC;AEE].');
end

