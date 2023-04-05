function [ti] = Find_Time_Match(te,Time_Matrix_Int)
% Time_Matrix_Int = [1 2 3 4;5 6 7 8;9 10 11 12];
% te = 6.6;
comp = double(te>Time_Matrix_Int);
temp = comp';
x = find(temp(:) == 0,1);
c = mod(x,size(Time_Matrix_Int,2)); %% Column
r = ((x-c)/size(Time_Matrix_Int,2))+1; %%Row
if c==0
    c=size(Time_Matrix_Int,2);
    r = r-1;
end
Post = Time_Matrix_Int(r,c);
if(c==1)
    if(r==1)
        Pre = Post; 
    else
        Pre = Time_Matrix_Int(r-1,size(Time_Matrix_Int,2));
    end
else
    Pre = Time_Matrix_Int(r,c-1);
end
if (abs(te-Post)>abs(te-Pre))
    ti = Pre;
else
    ti = Post;
end