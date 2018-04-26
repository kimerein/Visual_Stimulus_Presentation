function v = SeparateN(rrr)
% SeparateN takes a N-d matrix and separates it into N vectors
%
% v = SeparateN(rrr) takes a N-dim matrix with dimensions n1 x n2 x ... x nN and
% returns N vectors v{i}, i = 1:N, each with dimension ni x 1
%
% Note that each vector will have norm 1, so the amplitudes are
% meaningless.
%
% If N is 2, you are probably better off calling MakeSeparable so you also
% get a scaling factor, residuals, etc.
%
% If some of these dimensions are separable, then the relevant vectors are
% meaningful.
%
% For example, if the data are nstim x nt x nsites, and the different
% sites have different tuning, it does not make sense to look at v{1} and
% v{3}. But if all the sites and all the stimuli give responses with the
% same time course, then v{2} is meaningful.
%
% By the way, it probably makes most sense to remove the grand mean before
% calling this function.
%
% See also: MakeSeparable

% 2007-08 Matteo Carandini

%% Set up

nd = ndims(rrr);
nn = size(rrr);
n = numel(rrr);

if nd<2
    error('This function requires at least 2-dim input');
end

%% Do the work

v = {};
for id = 1:nd
    sss = shiftdim(rrr,id-1); % now the first dimension is the id
    rr = reshape(sss,nn(id),n/nn(id));
    [foo,v{id}] = MakeSeparable(rr);
end
clear foo

%% Graphics

figure; ax = [];
for id = 1:nd
    ax(id) = subplot(nd,1,id);
    plot(v{id});
end

%% Done

return







%% code to test the function (2 dims)

n1 = 10;
n2 = 20;

r = {};
r{1} = gaussian(2,1:n1);
r{2} = gaussian(2,(1:n2)-n2/2);

nd = length(r);

figure; ax = [];
for id = 1:nd
    ax(id) = subplot(nd,1,id);
    plot(r{id});
end

rrr = zeros(n1,n2);
for i1 = 1:n1
    for i2 = 1:n2
        rrr(i1,i2) = r{1}(i1)*r{2}(i2);
    end
end

v = SeparateN(rrr);

%% code to test the function (3 dims)

n1 = 10;
n2 = 20;
n3 = 30;

r = {};
r{1} = gaussian(2,1:n1);
r{2} = gaussian(2,(1:n2)-n2/2);
r{3} = sawtooth(1:n3);

nd = length(r);

figure; ax = [];
for id = 1:nd
    ax(id) = subplot(4,1,id);
    plot(r{id});
end

rrr = zeros(n1,n2,n3);
for i1 = 1:n1
    for i2 = 1:n2
        for i3 = 1:n3
            rrr(i1,i2,i3) = r{1}(i1)*r{2}(i2)*r{3}(i3);
        end
    end
end

rrr = rrr + + normrnd(0,.2,[n1 n2 n3]);

v = SeparateN(rrr);

%% code to test the function (4 dims)

n1 = 10;
n2 = 20;
n3 = 30;
n4 = 40;

r = {};
r{1} = gaussian(2,1:n1);
r{2} = gaussian(2,(1:n2)-n2/2);
r{3} = sawtooth(1:n3);
r{4} = square(1:n4);

nd = length(r);

figure; ax = [];
for id = 1:nd
    ax(id) = subplot(4,1,id);
    plot(r{id});
end

rrr = zeros(n1,n2,n3,n4);
for i1 = 1:n1
    for i2 = 1:n2
        for i3 = 1:n3
            for i4 = 1:n4
                rrr(i1,i2,i3,i4) = r{1}(i1)*r{2}(i2)*r{3}(i3)*r{4}(i4);
            end
        end
    end
end

rrr = rrr + + normrnd(0,.2,[n1 n2 n3 n4]);

v = SeparateN(rrr);

