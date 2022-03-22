generate_data( @diffusion1, 'diffusion1', 'lev', [2], 'deg', 2);
generate_data( @diffusion1, 'diffusion1', 'lev', [4], 'deg', 4);
generate_data( @diffusion1, 'diffusion1', 'lev', [5], 'deg', 5);

generate_data( @diffusion2, 'diffusion2', 'lev', [3,3], 'deg', 3);

function [ bcL, bcR ] = generate_data(pde_handle, pde_name, varargin)
data_dir = strcat("generated-inputs/boundary_conditions/");
out_format = strcat(data_dir, "compute_%s_%s_%dd_%dl_%dt_%ddim_%dp.dat");

root = get_root_folder();
[stat,msg] = mkdir ([root,'/gold/',char(data_dir)]);

%runtime_defaults
opts = OPTS(varargin);
opts.quiet = 1;
opts.use_oldcoeffmat = 0;
opts.use_oldhash = 0;
opts.max_lev_coeffs = 1;
%pde = check_pde(pde,opts);
pde = pde_handle( opts );

[elements, elements_idx]    = hash_table_sparse_nD (opts.lev, ...
                                                    opts.max_lev, opts.grid_type);
hash_table.elements         = elements;
hash_table.elements_idx     = elements_idx;

time = 0;

%%
% pde shortcuts

dims = pde.dimensions;
terms = pde.terms;

%%
% dim shortcuts

deg = opts.deg;
level = pde.dimensions{1}.lev;


num_terms = numel(pde.terms);
num_dims = numel(dims);

for tt = 1:num_terms % Construct a BC object for each term

    term_nD = terms{tt};

    for d1 = 1:num_dims

        dim = dims{d1};
        term_1D = term_nD.terms_1D{d1};
        %         type = term_1D_.type;

        xMin = dim.min;
        xMax = dim.max;

        lev = dim.lev;

        for p=1:numel(term_1D.pterms)
            this_type = term_1D.pterms{p}.type;
            this_g = term_1D.pterms{p}.g;
            this_dv = term_1D.pterms{p}.dV;

            if strcmp(this_type,'grad') || strcmp(this_type,'div')

                this_BCL = term_1D.pterms{p}.BCL;
                this_BCR = term_1D.pterms{p}.BCR;
                BCL_fList = term_1D.pterms{p}.BCL_fList;
                BCR_fList = term_1D.pterms{p}.BCR_fList;

                if strcmp(this_BCL,'D') % Left side

                    bcL = compute_boundary_condition(pde,this_g,this_dv,time,lev,deg,xMin,xMax,...
                                                     BCL_fList{d1},'L');
                    write_octave_like_output(sprintf(out_format, pde_name, "bcL", deg, level, ...
                                                     tt-1, d1-1, p-1), full(bcL)); 
                end

                if strcmp(this_BCR,'D') % Right side

                    bcR = compute_boundary_condition(pde,this_g,this_dv,time,lev,deg,xMin,xMax, ...
                                                     BCR_fList{d1},'R');

                    write_octave_like_output(sprintf(out_format, pde_name, "bcR", deg, level, ...
                                                     tt-1, d1-1, p-1), full(bcR)); 
                end

            end

        end

    end % loop over dim1

end % loop over terms


end
