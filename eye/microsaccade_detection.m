function ms = microsaccade_detection(eye_x, eye_y, samplingRate, algorithm, fig)
%% 
% detect mscade events in the given x,y eye-positions
% during fixation by using a specified algorithm
%
% INPUT: eye_x ... vector of horizontal eye position in one trial of fixation
%        eye_y ... vector of vertical eye position in one trial of fixation
%        samplingRate ... sampling rate in the measurement system
%        algorithm ... 'Nienborg'; Nienborg & Cumming (2006) ' http://www.jneurosci.org/content/26/37/9567.long' 
%                  ... 'Engert'; Engbert & Kliegl (2003) 'http://www.sciencedirect.com/science/article/pii/S0042698903000841'
%                  ... 'Joshi'; Joshi et al (2016) 'http://www.sciencedirect.com/science/article/pii/S089662731501034X'
%                  ... 'hafed'; Hafed et al (2011) 'http://www.jneurosci.org/content/31/43/15219.long'
%        fig ... 1; validation plot
% OUTPUT: ms ... matlab structure file storing relevant infos
%
% NOTE: As a detection algorithm 'Engert' is recommended, as this is the
% one widely used.
%
% EXAMPLE: ms = mscade(eye_x, eye_y, 500, 'Nienborg', 1);
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++

%%
% deal with inputs
if nargin < 2 
    error('At least eye-positions for x and y must be given as the first and second input arguments.')
end
if nargin < 3; samplingRate = 500; end
if nargin < 4; algorithm = 'Nienborg'; end
if nargin < 5; fig = 1; end

%%
% match the length of eye_x and eye_y, if different
if length(eye_x) ~= length(eye_y)
  disp('eye-position x and y have different length.')
  % force eye_x and eye_y to be the same length
    while abs(length(eye_x)-length(eye_y))>0
        if length(eye_x) > length(eye_y)
            eye_x(end) = [];
        else
            eye_y(end) = [];
        end
    end
    disp('Length of eye-position x and y were matched.')
end
len_eye = length(eye_x);

%%
% define a detection algorithm
switch lower(algorithm)
    case  'nienborg' % Nienborg & Cumming (2006)
        ms.algorithm.name = 'Nienborg & Cumming (2006)';
        ms.algorithm.velocity = 12; % deg/s
        ms.algorithm.smooth = 0;
        ms.algorithm.acceleration = nan;
        ms.algorithm.minduration = nan;

      case 'engbert' % Engbert & Kliegl (2003); Engbert & Mergenthaler (2006)
        % transform eye-positions into velocities
        vel_x = zeros(1,len_eye);
        vel_y = zeros(1,len_eye);
        dt = 1/samplingRate;
        for i = 3:len_eye-2
              vel_x(i) = (eye_x(i+2)+eye_x(i+1)-eye_x(i-1)-eye_x(i-2))...
                    /(6*dt);
              vel_y(i) = (eye_y(i+2)+eye_y(i+1)-eye_y(i-1)-eye_y(i-2))...
                    /(6*dt);
        end

        % SD of the velocity distribution as the detection threshold
        sigma = zeros(2,1);
        sigma(1) = median(vel_x.^2) - (median(vel_x))^2;
        sigma(2) = median(vel_y.^2) - (median(vel_y))^2;
        gamma = 6; % default used in the original paper
        thre_x = gamma*sigma(1);
        thre_y = gamma*sigma(2);
        
        % store algorithm info
        ms.algorithm = 'Engbert & Kliegl (2003)';
        ms.algorithm.velocity = [thre_x, thre_y]; % deg/s
        ms.algorithm.smooth = 1;
        ms.algorithm.acceleration = nan;
        ms.algorithm.minduration = 6; % ms
            
    case 'joshi' % Joshi et al., 2016
        ms.algorithm = 'Joshi et al (2016)';
        ms.algorithm.velocity = 15; % deg/s
        ms.algorithm.smooth = 0;
        ms.algorithm.acceleration = nan;
        ms.algorithm.minduration = 6; % ms
        
   case 'hafed' % Hafed et al., 2011        
        ms.algorithm = 'Hafed et al (2011)';
        ms.algorithm.velocity = 8; % deg/s
        ms.algorithm.smooth = 0;
        ms.algorithm.acceleration = 550;
        ms.algorithm.minduration = nan; % ms
end

%%
% detection ms by velocity
event = zeros(1,len_eye);
if ms.algorithm.smooth == 0
    % compute instantaneous eye velocity
    vel_x = [0 diff(eye_x)*samplingRate];
    vel_y = [0 diff(eye_y)*samplingRate];
    vel = sqrt(vel_x.^2 + vel_y.^2);
    event(vel > ms.algorithm.velocity) = 1;      
elseif ms.algorithm.smooth == 1
    % detect ms events           
    event((vel_x/ms.algorithm.velocity(1)).^2 + ...
        (vel_y/ms.algorithm.velocity(2)).^2 > 1) = 1;
end

% detection ms by acceleration
if ~isnan(ms.algorithm.acceleration)
    vel = sqrt(vel_x.^2 + vel_y.^2);
    acc = [0 diff(vel)];
    event(acc > ms.algorithm.acceleration) = 1;
end

% detection ms by minimum duration
if ~isnan(ms.algorithm.minduration)
    minsamples = 1;
else
    minsamples = mindur.minduration*dt*1000;
end

% apply minimum duration criteria
pos = 1;
sacc_start = [];
sacc_end = [];
while pos < len_eye
   idx = find(event(pos:end)==1, 1, 'first');
   s = 1;
   stop = 0;
   while stop==0
       if event(idx+s)==1
           s = s + 1;
       else
           s = s - 1;
           stop = 1;
       end
   end
   if s < minsamples
       event(idx:idx+s) = 0;
   else
       sacc_start = [sacc_start, idx];
       sacc_end = [sacc_end, idx+s];
   end
   pos = idx + s + 1;
end

% store eye-data into a structure
if isempty(sacc_start)
    ms.counts = 0;
else
    ms.counts = length(sacc_start);
end
if ms.counts > 0
    ms.amp = zeros(1, ms.counts);
    ms.peakv = zeros(1, ms.counts);
    ms.duration = zeros(1, ms.counts);
    ms.angle = zeros(1, ms.counts);
    for i = 1:length(sacc_start)
        ms.amp(i) = sqrt((eye_x(sacc_end(i)) - eye_x(sacc_start(i))).^2 ...
            + (eye_y(sacc_end(i)) - eye_y(sacc_start(i))).^2);
        ms.peakv(i) = max(sqrt(diff(eye_x(sacc_start(i):sacc_end(i))).^2 ...
            + diff(eye_y(sacc_start(i):sacc_end(i))).^2))*samplingRate;
        ms.duration(i) = (sacc_end(i) - sacc_start(i))/samplingRate;
        ms.angle(i) = atan2(eye_y(sacc_end(i)) - eye_y(sacc_start(i)), ...
            eye_x(sacc_end(i)) - eye_x(sacc_start(i)))*180/pi;
    end
else
    ms.amp = nan;
    ms.peakv = nan;
    ms.duration = nan;
    ms.angle = nan;
end
ms.velocity_x = vel_x;
ms.velocity_y = vel_y;
ms.event = event;

% validate with plot
if fig==1
      figure;
      % time domain
      subplot(4,6,1:6)
      time = [1:length(eye_x)]/samplingRate;
      plot(time, eye_x, '-','Color',[0.5 0.5 0.5])
      ylabel('x (deg)')
      set(gca,'box','off'); set(gca,'TickDir','out');  
      subplot(4,6,7:12)
      plot(time, eye_y, '-','Color',[0.5 0.5 0.5])
      ylabel('y (deg)')
      xlabel('time (sec)')
      set(gca,'box','off'); set(gca,'TickDir','out');  

      % position space
      subplot(4,6,13:14)
      plot(eye_x,eye_y,'-','Color',[0.5 0.5 0.5])
      xlabel('horizontal position (deg)')
      ylabel('vertical position (deg)')
      set(gca,'box','off'); set(gca,'TickDir','out');  

      % velocity space
      subplot(4,6,17:18)
      if ms.algorithm.smooth == 1
        thre = sqrt(thre_x^2 + thre_y^2);
      end
      th = 0:pi/50:2*pi;
      xunit = thre * cos(th);
      yunit = thre * sin(th);
      plot(xunit, yunit,'--k','lineWidth',0.5);
      hold on;
      plot(vel_x, vel_y, '-','Color',[0.5 0.5 0.5])
      xlabel('horizontal velocity (deg/s)')
      ylabel('vertical velocity (deg/s)')
      set(gca,'box','off'); set(gca,'TickDir','out');  

      % draw saccadic traces
      map = lines(length(sacc_start));
      for i = 1:length(sacc_start)
          k = 0;
          while sacc_start(i) + k < sacc_end(i)
              subplot(4,6,1:6)
              hold on;
              plot(time(sacc_start(i)+k:sacc_start(i)+k+1), eye_x(sacc_start(i)+k:sacc_start(i)+k+1),...
                  '-','Color',map(i,:))
              subplot(4,6,7:12)
              hold on;
              plot(time(sacc_start(i)+k:sacc_start(i)+k+1), eye_y(sacc_start(i)+k:sacc_start(i)+k+1),...
                  '-','Color',map(i,:))
              subplot(4,6,13:15)
              hold on;
              plot([eye_x(sacc_start(i)+k) eye_x(sacc_start(i)+k+1)], [eye_y(sacc_start(i)+k) eye_y(sacc_start(i)+k+1)],...
                      '-','Color',map(i,:),'lineWidth',1.5)
              subplot(4,6,16:18)
              hold on;
              plot([vel_x(sacc_start(i)+k) vel_x(sacc_start(i)+k+1)], [vel_y(sacc_start(i)+k) vel_y(sacc_start(i)+k+1)],...
                      '-','Color',map(i,:),'lineWidth',1.5)
              k = k + 1;
          end
      end

      subplot(4,6,[19 20])
      plot(ms.amp, ms.peakv,'ok')
      set(gca, 'XScale', 'log')
      set(gca, 'YScale', 'log')
      xlabel('amplitude (deg)')
      ylabel('peak velocity (deg/s)')
      set(gca,'box','off'); set(gca,'TickDir','out');  

      subplot(4,6,[21 22])
      polarhistogram(ms.angle*pi/180)
      set(gca,'box','off'); set(gca,'TickDir','out'); 
      
      subplot(4,6,[23 24])
      histogram(1000*ms.duration')
      xlabel('duration (ms)')
      set(gca,'box','off'); set(gca,'TickDir','out'); 
end
