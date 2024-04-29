classdef RoomWithAAES
   properties
      room_dims {mustBeNumeric};
      alpha_set {mustBeNumeric};
      src_positions {mustBeNumeric};
      rec_positions {mustBeNumeric};
      ls_positions {mustBeNumeric};
      mic_positions {mustBeNumeric};
      sample_rate {mustBeNumeric};
      bit_depth {mustBeNumeric};
   end
   methods
       % Constructor
       function obj = RoomWithAAES(room_dims, alpha_set, src_positions, rec_positions, ls_positions, mic_positions, sample_rate, bit_depth)
           obj.room_dims = room_dims;
           obj.alpha_set = alpha_set;
           obj.src_positions = src_positions;
           obj.rec_positions = rec_positions;
           obj.ls_positions = ls_positions;
           obj.mic_positions = mic_positions;
           obj.sample_rate = sample_rate;
           obj.bit_depth = bit_depth;
       end
       
       % Generates the room impulse responses for each
       % loudspeaker-microphone pair in the room
       function GenerateSystemIRs(obj, output_parent_dir)
           output_dir = output_parent_dir + "AAES IRs Ch["+size(obj.mic_positions, 1)+"x"+size(obj.ls_positions, 1)+"] Room"+mat2str(obj.room_dims)+" AlphaSet["+obj.alpha_set+"]/";
           mkdir(output_dir);

           alphas = readmatrix("Example Absorption Coefficients/alpha_set_" + obj.alpha_set + ".dat");

           should_normalise = true; % This will batch normalise all outputs to 0 dBFS, preserving level relationships

           GenerateAKToolsRIRs(obj.room_dims, alphas, obj.src_positions, obj.rec_positions, obj.ls_positions, obj.mic_positions, obj.sample_rate, output_dir, obj.bit_depth, should_normalise);
      end
   end
end