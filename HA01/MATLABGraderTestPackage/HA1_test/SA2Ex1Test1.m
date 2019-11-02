P_D = rand;
lambda_c = rand;
range_c = [-rand rand;-rand rand];

% Run learner solution.
sensor_model = modelgen.sensormodel(P_D,lambda_c,range_c);

% Run reference solution.
sensor_model_ref = reference.modelgen.sensormodel(P_D,lambda_c,range_c);

% Compare.
% assessVariableEqual('sensor_model', sensor_model_ref);
tol = 1e-4;
assert(abs(sensor_model.pdf_c-sensor_model_ref.pdf_c) < tol,'The value of clutter pdf 2D case: Expected %f, got %f.', sensor_model_ref.pdf_c, sensor_model.pdf_c);
assert(abs(sensor_model.intensity_c-sensor_model_ref.intensity_c) < tol,'Clutter intensity 2D case: Expected %f, got %f.', sensor_model_ref.intensity_c, sensor_model.intensity_c);

range_c = [-rand rand];

% Run learner solution.
sensor_model = modelgen.sensormodel(P_D,lambda_c,range_c);

% Run reference solution.
sensor_model_ref = reference.modelgen.sensormodel(P_D,lambda_c,range_c);

% Compare.
% assessVariableEqual('sensor_model', sensor_model_ref);
tol = 1e-4;
assert(abs(sensor_model.pdf_c-sensor_model_ref.pdf_c) < tol,'The value of clutter pdf 1D case: Expected %f, got %f.', sensor_model_ref.pdf_c, sensor_model.pdf_c);
assert(abs(sensor_model.intensity_c-sensor_model_ref.intensity_c) < tol,'Clutter intensity 1D case: Expected %f, got %f.', sensor_model_ref.intensity_c, sensor_model.intensity_c);

