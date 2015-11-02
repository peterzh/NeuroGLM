%% Simulate from a given model
expRef = '2015-10-30_1_Hopkins'; %Used for getting input contrasts from only
g = GLM(expRef).setModel('C^N');

%Set parameter values & simulate responses (and insert into data fields)
sim_params = [0 10 1 0 1 10 0.5];
sim_phat = g.calculatePhat(sim_params,g.data.contrast_cond);

choice = @(phat)(sum(rand>[0,cumsum(phat(1:2))]));
sim_choices = [];
for i = 1:length(sim_phat)
    g.data.response(i,1) = choice(sim_phat(i,:));
end

g.data.repeatNum = ones(length(sim_phat),1);
g.expRef = '2015-11-02_1_SIMULATION';

%% Run each model crossvalidated fitting and save results
saveDir = '\\basket.cortexlab.net\homes\peterzh\NeuroGLM\ModelFiles';
models = {'Offset','ifC','C^N','C^N-subset','C^NL^NR','C^NL^NR-subset','C50','C50-subset','Supersaturation-subset'};

for m = 1:length(models)
    g = g.setModel(models{m});
    g = g.fit('crossval');
    save(fullfile(saveDir,[g.expRef '_' g.modelString '_crossvalidated.mat']),'g');
    disp('done');
end

