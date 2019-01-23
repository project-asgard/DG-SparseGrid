function Y = krond(nkron, Acell)
% Y = krond(nkron, Acell)
% n-dimensional version of kron
%
% Y = kron( Acell{1}, Acell{2}, ..., Acell{nkron} )
%
if (nkron == 1),
    Y = Acell{1};
    return;
end

if (nkron == 2),
    Y = kron( Acell{1}, Acell{2} );
    return;
end

Y = kron(krond( nkron-1, Acell), Acell{nkron});
end
