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
% EXAMPLE: ms = mscade(eye_x, eye_y, 500, 'Nienborg', 1);
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++

if nargin < 2 
    error('At least eye-positions for x and y must be given as the first and second input arguments.')
end
if nargin < 3; samplingRate = 500; end
if nargin < 4; algorithm = 'Nienborg'; end
if nargin < 5; fig = 1; end

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

switch lower(algorithm)
    case  'nienborg' % Nienborg & Cumming (2006)
        n = length(eye_x);
        % compute saccadic velocity
        vel_x = [0 diff(eye_x)*samplingRate];
        vel_y = [0 diff(eye_y)*samplingRate];
        vel = sqrt(vel_x.^2 + vel_y.^2);
        
        % start: 1, end: -1
        thre = 12; % degree / sec
        event = zeros(1,n);
        event(vel > thre) = 1;      
        diff_event = [0 diff(event)];
        sacc_start = find(diff_event==1);
        sacc_end = find(diff_event==-1) - 1;

        % equalize start and end
        delta = abs(length(sacc_start) - length(sacc_end));
        if length(sacc_start) > length(sacc_end)
                sacc_start(end-delta+1:end) = [];
        elseif length(sacc_start) < length(sacc_end)
                sacc_end(end-delta+1:end)  = [];
        end
        
        % remove saccade events with amplitude 0
        zeroamp = sacc_start==sacc_end;
        sacc_start(zeroamp) = [];
        sacc_end(zeroamp) = [];
        
        % equalize start and end
        delta = abs(length(sacc_start) - length(sacc_end));
        if length(sacc_start) > length(sacc_end)
                sacc_start(end-delta+1:end) = [];
        elseif length(sacc_start) < length(sacc_end)
                sacc_end(end-delta+1:end)  = [];
        end
               
        % store basic parameters
        ms.algorithm = 'Nienborg & Cumming (2006)';
        ms.threshold.velocity = thre;

      case 'engbert' % Engbert & Kliegl (2003); Engbert & Mergenthaler (2006)
            % transform eye-positions into velocities
            len_eye = length(eye_x);
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
                      
            % detect mscade events           
            event = zeros(1,len_eye);
            event((vel_x/thre_x).^2 + (vel_y/thre_y).^2 > 1) = 1;
            diff_event = [0 diff(event)];

            % amplitude, peak velocity, duration
            sacc_start = find(diff_event==1) ;
            sacc_end = find(diff_event==-1) - 1; 
            
            % remove saccade events with amplitude 0
            zeroamp = find(sacc_start==sacc_end);
            sacc_start(zeroamp) = [];
            sacc_end(zeroamp) = [];
                        
            % equalize start and end
            delta = abs(length(sacc_start) - length(sacc_end));
            if length(sacc_start) > length(sacc_end)
                sacc_start(end-delta+1:end) = [];
            elseif length(sacc_start) < length(sacc_end)
                sacc_end(end-delta+1:end)  = [];
            end
        
            % store basic info
            ms.algorithm = 'Engbert & Kliegl (2003)';
            ms.threshold.velocity = [thre_x, thre_y];
            
    case 'joshi' % Joshi et al., 2016
        n = length(eye_x);
        % compute saccadic velocity
        vel_x = [0 diff(eye_x)*samplingRate];
        vel_y = [0 diff(eye_y)*samplingRate];
        vel = sqrt(vel_x.^2 + vel_y.^2);
        
        % start: 1, end: -1
        thre = 15; % degree / sec
        event = zeros(1,n);
        event(vel >= thre) = 1;      
        diff_event = [0 diff(event)];
        sacc_start = find(diff_event==1);
        sacc_end = find(diff_event==-1) - 1;
        
       % remove saccade events with amplitude 0
        zeroamp = find(sacc_start==sacc_end);
        sacc_start(zeroamp) = [];
        sacc_end(zeroamp) = [];
        
        % equalize start and end
        delta = abs(length(sacc_start) - length(sacc_end));
        if length(sacc_start) > length(sacc_end)
                sacc_start(end-delta+1:end) = [];
        elseif length(sacc_start) < length(sacc_end)
                sacc_end(end-delta+1:end)  = [];
        end
        
        % criteria of duration
        for l = 1:length(sacc_start)
            dur = (sacc_end(l) - sacc_start(l) + 1)/samplingRate;
            if dur < 0.006
                event(sacc_start(l):sacc_end(l)) = 0;
            end
        end

        % redefine event
        diff_event = [0 diff(event)];        
        
        % start and end of saccade
        sacc_start = find(diff_event==1);
        sacc_end = find(diff_event==-1) -1;
        
        % equalize start and end again
        delta = abs(length(sacc_start) - length(sacc_end));
        if length(sacc_start) > length(sacc_end)
                sacc_start(end-delta+1:end) = [];
        elseif length(sacc_start) < length(sacc_end)
                sacc_end(end-delta+1:end)  = [];
        end
        
        % store basic info
        ms.algorithm = 'Joshi et al (2016)';
        ms.threshold.velocity = thre;
        ms.threshold.duration = 0.006;
        
   case 'hafed' % Hafed et al., 2011
        n = length(eye_x);
        % compute saccadic velocity
        vel_x = [0 diff(eye_x)*samplingRate];
        vel_y = [0 diff(eye_y)*samplingRate];
        vel = sqrt(vel_x.^2 + vel_y.^2);
        
        % compute acceleration
        acc_x = [0 diff(vel_x)];
        acc_y = [0 diff(vel_y)];
        acc = sqrt(acc_x.^2 + acc_y.^2);
                
        % start: 1, end: -1
        thre = 8; % degree / sec
        thre2 = 550; % degree /sec^2
        event = zeros(1,n);
        event(vel >= thre & acc > thre2) = 1;      
        diff_event = [0 diff(event)];
        sacc_start = find(diff_event==1);
        sacc_end = find(diff_event==-1) - 1;
        
       % remove saccade events with amplitude 0
        zeroamp = find(sacc_start==sacc_end);
        sacc_start(zeroamp) = [];
        sacc_end(zeroamp) = [];
        
        % equalize start and end
        delta = abs(length(sacc_start) - length(sacc_end));
        if length(sacc_start) > length(sacc_end)
                sacc_start(end-delta+1:end) = [];
        elseif length(sacc_start) < length(sacc_end)
                sacc_end(end-delta+1:end)  = [];
        end
                
        % store basic info
        ms.algorithm = 'Hafed et al (2011)';
        ms.threshold.velocity = thre;
        ms.threshold.acceleration = thre2;
end

% store eye-data into a structure
ms.counts = length(sacc_start);
ms.amp = zeros(1, ms.counts);
ms.peakv = zeros(1, ms.counts);
ms.duration = zeros(1, ms.counts);
ms.angle = zeros(1, ms.counts);
for i = 1:length(sacc_start)
    k = 0;
    temp_amp = [];
    temp_pv = [];
    angle = [];
    while sacc_start(i) + k < sacc_end(i)
        vx = eye_x(sacc_start(i)+k+1) - eye_x(sacc_start(i)+k);
        vy = eye_y(sacc_start(i)+k+1) - eye_y(sacc_start(i)+k);
        angle = [angle -atan2(vy, vx)*180/pi]; % note the negative y
        temp_amp = [temp_amp sqrt(vx^2 + vy^2)];
        temp_pv = [temp_pv max(sqrt(vel_x(sacc_start(i)+k:sacc_start(i)+k+1).^2 + vel_y(sacc_start(i)+k:sacc_start(i)+k+1).^2))];
        k = k + 1;
    end
    ms.amp(i) = sum(temp_amp);
    ms.peakv(i) = max(temp_pv);
    ms.duration(i) = k/samplingRate;
    ms.angle(i) = median(angle);
end
            
ms.eye_x = eye_x;
ms.eye_y = eye_y;
ms.velocity_x = vel_x;
ms.velocity_y = vel_y;
ms.event = diff_event;

% validate with plot
if fig==1
      figure;
      % time domain
      subplot(4,6,1:6)
      time = [1:length(eye_x)]/samplingRate;
      plot(time, eye_x, '-','Color',[0.5 0.5 0.5])
      ylabel('x (deg)')
      set(gca,'box','off'); set(gca,'TickDir','out'); axis square 
      subplot(4,6,7:12)
      plot(time, eye_y, '-','Color',[0.5 0.5 0.5])
      ylabel('y (deg)')
      xlabel('time (sec)')
      set(gca,'box','off'); set(gca,'TickDir','out'); axis square 

      % position space
      subplot(4,6,13:15)
      plot(eye_x,eye_y,'-','Color',[0.5 0.5 0.5])
      xlabel('horizontal position (deg)')
      ylabel('vertical position (deg)')
      set(gca,'box','off'); set(gca,'TickDir','out'); axis square 

      % velocity space
      subplot(4,6,16:18)
      switch lower(algorithm)
          case 'engbert'
              plot([-thre_x thre_x],[-thre_y -thre_y],'--k','lineWidth',0.5)
              hold on;
              plot([-thre_x thre_x],[thre_y thre_y],'--k','lineWidth',0.5)
              hold on;
              plot([-thre_x -thre_x],[-thre_y thre_y],'--k','lineWidth',0.5)
              hold on;
              plot([thre_x thre_x],[-thre_y thre_y],'--k','lineWidth',0.5)
          otherwise
              th = 0:pi/50:2*pi;
              xunit = thre * cos(th);
              yunit = thre * sin(th);
              plot(xunit, yunit,'--k','lineWidth',0.5);
      end
      hold on:
      plot(vel_x, vel_y, '-','Color',[0.5 0.5 0.5])
      xlabel('horizontal velocity (deg/s)')
      ylabel('vertical velocity (deg/s)')
      set(gca,'box','off'); set(gca,'TickDir','out'); axis square 

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
      plot(ms.amp, ms.peakv,'ok'),lsline
      xlabel('amplitude (deg)')
      ylabel('peak velocity (deg/s)')
      set(gca,'box','off'); set(gca,'TickDir','out'); axis square 

      subplot(4,6,[21 22])
      polarhistogram(ms.angle*pi/180)
      set(gca,'box','off'); set(gca,'TickDir','out'); axis square
      
      subplot(4,6,[23 24])
      histogram(1000*ms.duration')
      xlabel('duration (ms)')
      set(gca,'box','off'); set(gca,'TickDir','out'); axis square
end
