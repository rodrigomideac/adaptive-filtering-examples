function assertDimensions(A,B)
    assert(isequal(size(A), size(B)) || (isvector(A) && isvector(B) && numel(A) == numel(B)),'Dimensions dont match')
end
