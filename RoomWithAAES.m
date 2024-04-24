classdef RoomWithAAES
   properties
      room_dims {mustBeNumeric};
      alphas {mustBeNumeric};
      src_positions {mustBeNumeric};
      rec_positions {mustBeNumeric};
      ls_positions {mustBeNumeric};
      mic_positions {mustBeNumeric};
      sample_rate {mustBeNumeric};
      bit_depth {mustBeNumeric};
   end
   methods
       % Constructor
       function obj = RoomWithAAES(room_dims, alphas, src_positions, rec_positions, ls_positions, mic_positions, sample_rate, bit_depth)
           obj.room_dims = room_dims;
           obj.alphas = alphas;
           obj.src_positions = src_positions;
           obj.rec_positions = rec_positions;
           obj.ls_positions = ls_positions;
           obj.mic_positions = mic_positions;
           obj.sample_rate = sample_rate;
           obj.bit_depth = bit_depth;
       end
       
       % Generates the room impulse responses for each
       % loudspeaker-microphone pair, for each room, for each set of
       % absorption coefficients
       function GenerateRIRs(obj)
           for room_num = 1:size(obj.room_dims, 1)
               for alpha_set = 1:size(obj.alphas, 1)
                   output_dir = "Automated RIRs/AAES IRs Ch["+size(obj.mic_positions, 1)+"x"+size(obj.ls_positions, 1)+"] Room"+mat2str(obj.room_dims(room_num, :))+" AlphaSet["+alpha_set+"]/";
                   mkdir(output_dir);

                   should_normalise = true; % This will batch normalise all outputs to 0 dBFS, preserving level relationships
        
                   GenerateRIRs(obj.room_dims, obj.alphas, obj.src_positions, obj.rec_positions, obj.ls_positions, obj.mic_positions, obj.sample_rate, output_dir, obj.bit_depth, should_normalise);
               end
           end
      end
   end
end