%% Peak Detector - basic (finding all peaks)

%[y,Fs] = audioread('DSP_Project_Test.wav');
[y,Fs] = audioread('HappyBirthday_test_sharp.wav');
dim = size(y);

samples = [];

if(dim(2)==2)
    samples = y(:,1);
else
    samples = y;
end


%sample = y(:); % 1D matrix
indices = [];

Fs

for i = 1:size(samples)
   indices(i) = i;
end

mpprom = 0.4 * max(samples);
[peaks, locs] = findpeaks(samples, indices, 'MinPeakHeight', mpprom);

findpeaks(samples, indices, 'MinPeakHeight', mpprom)
xlabel('Sample')
ylabel('Index?')
title('Find All Peaks')

peaks
locs = locs';
locs

distances_btw = []; % stores the distances between peaks

for i = 1:(size(locs)-1)
    distances_btw(i, 1) = locs(i+1) - locs(i);
end

distances_btw

long_distances = [];
locs_filtered = [];

j = 1;
for i = 1:size(distances_btw)
    if distances_btw(i) > 2000 % arbitrary
        long_distances(j) = distances_btw(i);
        locs_filtered(j) = locs(i);
        j=j+1;
    end    
end
locs_filtered(j) = locs(i+1);


bpm_distance = mean(long_distances);

averages = [0,0,0];
locs_filtered = locs_filtered';
% find silence value
for i=1:size(locs_filtered)-1
    middle = locs_filtered(i) + (long_distances(i) ./ 2);
    low = middle - long_distances(i) ./ 4;
    high = middle + long_distances(i) ./ 4;
    avg = mean(abs(samples(low:high)));
    averages(i) = avg;
end

silence = mean(averages);

eighth_beat_indices = [];
eighth_beat_values = [];
j = 1;
i = locs_filtered(size(locs_filtered));

disp('-------------');
i = (i(1,1)) + bpm_distance;
y_size = size(samples(:,1));
y_size = y_size(1,1);
disp(i);
disp(y_size);
disp('-------------');

% i holds the sample number of the snap location

while i < y_size
    %disp(i);
    eighth_beat_indices(j) = indices(int32(i));
    eighth_beat_values(j) = samples(int32(i));
    i = i + bpm_distance ./ 2;
    j = j+1;
end

eighth_beat_indices = eighth_beat_indices';
eighth_beat_values = eighth_beat_values';

i = 1;
j = 1;

binary_vector = [size(eighth_beat_values)];

note_durations = zeros([], 2);
k = 1;
note_playing = 0;
rest_playing = 0;

%disp('here')

num_eighth_beats = size(eighth_beat_values);
num_eighth_beats = num_eighth_beats(1,1);

beat_indexer = locs_filtered(end) + bpm_distance;
last_max_index = locs_filtered(end);
bpm_adjusted = bpm_distance / 2;

mpheight = 0.6 * max(samples(beat_indexer:end));

while beat_indexer < y_size    % Look over range of 1/16 note around where note SHOULD be
    disp('----------');
   
    disp(beat_indexer);
    lowerbound = int32(beat_indexer - (bpm_distance ./4));
    upperbound = int32(beat_indexer + (bpm_distance ./4));
    if(upperbound > y_size)
        upperbound = y_size;
    end
    
    % Will substitute median for mean 
    
    sorted = sort(abs(samples(lowerbound:upperbound)));
    sortedSize = size(sorted);
    median = sorted(int32(sortedSize(1)/2));
    
    % Find if there is a peak
    peaks2=[]; locs2=[];
    [peaks2, locs2] = findpeaks(samples(lowerbound:upperbound), indices(lowerbound:upperbound), 'MinPeakHeight', mpheight);
    %figure
    %findpeaks(samples(lowerbound:upperbound), indices(lowerbound:upperbound), 'MinPeakHeight', mpheight)
    
    [max_value, max_value_loc] = max(samples(lowerbound:upperbound));
    
    fprintf('    lowerbound = %i', lowerbound);
    fprintf('    upperbound = %i', upperbound);
    fprintf('    median     = %i', median);
    fprintf('    max_value  = %i', max_value);
    fprintf('    max_index  = %i', max_value_loc + lowerbound);
    fprintf('');
    %disp('median / max');
    %median
    %max_value
    
    if size(peaks2,1) > 0 | (k > 1 & median >= silence * 1.5 & note_playing == 0)
        % new note
        disp('new note');
        if(note_playing == 1)
            note_durations(k-1, 2) = beat_indexer;
        end
        note_durations(k, 1) = beat_indexer;
        note_playing = 1;
        k = k+1;
        bpm_adjusted = ((max_value_loc+lowerbound) - last_max_index)/2;
    elseif median < silence * 1.5
        % rest
        disp('theres a rest');
        if note_playing == 1
            note_durations(k, 2) = beat_indexer;
            note_playing = 0;
        else
            ;;
        end
        k = k + 1;
    else
        disp('note continued');
        % continued note
    end
    
    %bpm_adjusted = (bpm_adjusted + bpm_distance/2) / 2;
    bpm_adjusted = bpm_distance/2;
    last_max_index = beat_indexer;
    beat_indexer = beat_indexer + bpm_adjusted;
    disp(beat_indexer);
end

disp('here')

note_durations