function result = binary_multiplication(a,b)
    c = a*b;
    q = size(c);
    rows = q(1);
    columns = q(2);
    
    for i =1:rows
        for j = 1:columns
            if mod(c(i,j),2)==0
                result(i,j) = 0;
            end
            if mod(c(i,j),2) > 0
                result(i,j) = 1;
            end
        end
    end

end