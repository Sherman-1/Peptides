    

def elongate_tm_segments(tm_indices : dict, protein_length : int, min_length=20, max_length=70):
    """
    This function takes a list of tuples containing the start and end indices of putative transmembrane (tm) segments
    Extracted for the same multiple-fragments transmembrane protein.
    For example, GPCR proteins have 7 transmembrane segments, they will end up in a list of 7 tuples.

    For each tm segment, the function will elongate the segment to a random size drawn from a given size distribution,
    given by the size_picker_v2 function.

    The function will elongate the segment only if the size of the segment is smaller than the size drawn from the distribution.
    The goal here is to "draw" from the parts of the sequence that are not transmembrane segments, and elongate the tm segments.
    One main goal is to avoid drawing twice from the same region to elongate two tm segments that are adjacent to each other.

    Input:

    tm_indices[chain_id] : list of tuples 
                # [ (12,26,15), (45, 60, 16), (80, 100, 21) ...]
                # [ (start, end, length), ... ]

    min_length : int
                # minimum length of the elongated segment

    max_length : int
                # maximum length of the elongated segment
    """

    for chain_id in tm_indices:

        

        ##### Treat first TM Segment separately ##### 


        desired_length = size_picker_v2(min_size=min_length, max_size=max_length)

        
        # First TM Segment
        start_current = tm_indices[chain_id][0][0]
        end_current = tm_indices[chain_id][0][1]
        length_current = tm_indices[chain_id][0][2]


        if desired_length > length_current:


            # Second TM Segment
            start_next = tm_indices[chain_id][1][0]

            elongation_left_to_do = desired_length - length_current


            downstream = random.randint(0, elongation_left_to_do)

            
            lefts = None

            # The new end of this tm segment should not exceed the start of the next tm segment
            if downstream + end_current > start_next:

                new_end_coordinates = start_next - 1

                lefts = downstream - ( start_next - end_current )



            else:

                new_end_coordinates = end_current + downstream


            upstream = elongation_left_to_do - downstream

            

            if lefts:

                upstream += lefts



            if start_current - upstream < 1:

                new_start_coordinates = 1



            else:

                new_start_coordinates = start_current - upstream

            tm_indices[chain_id][0] = (new_start_coordinates, new_end_coordinates, new_end_coordinates - new_start_coordinates)


        ##### Treat from the second TM Segment to the penultimate one ( n-1 ) #####

        for i in range(1, len(tm_indices[chain_id]) - 1):

            # Target size that the current tm should reach
            desired_length = size_picker_v2(min_size=min_length, max_size=max_length)

            # ith TM Segment
            start_current = tm_indices[chain_id][i][0]
            end_current = tm_indices[chain_id][i][1]
            length_current = tm_indices[chain_id][i][2]

            # check before anything else to save computation time
            if desired_length <= length_current:

                # If there is no elongation to do, we skip to the next segment
                # and the coordinates of the ith segment are not modified
                continue
            
            # (i+1)th TM Segment
            start_next = tm_indices[chain_id][i+1][0]


            # (i-1)th TM Segment
            end_previous = tm_indices[chain_id][i-1][1]
            
            # Compute the number of residues that are required to elongate the current segment
            elongation_left_to_do = desired_length - length_current


            # Randomly choose the number of residues to elongate downstream ( toward the C-terminal )
            downstream = random.randint(0, elongation_left_to_do)

            lefts = None

            # The new end of this tm segment should not exceed the start of the next tm segment
            if downstream + end_current > start_next:

                # Hence take everyting that is between the end of the current tm segment and the start of the next one
                new_end_coordinates = start_next - 1

                # What is " left " from downstream that could not be taken cause of the next tm ? 
                lefts = downstream - (start_next - end_current)

            else:

                new_end_coordinates = end_current + downstream

            ## If there is elongation that was not taken from downstream, add it to the upstream
            upstream = elongation_left_to_do - downstream
            if lefts:

                upstream += lefts


            # The new start of this tm segment should not be lower than the end of the previous tm segment
            if start_current - upstream < end_previous:

                new_start_coordinates = end_previous + 1 

            else:

                new_start_coordinates = start_current - upstream


            tm_indices[chain_id][i] = (new_start_coordinates, new_end_coordinates, new_end_coordinates - new_start_coordinates)

            


        ##### Treat the last TM Segment #####

        # Target size that the current tm should reach
        desired_length = size_picker_v2(min_size=min_length, max_size=max_length)


        # Last TM Segment
        start_current = tm_indices[chain_id][-1][0]
        end_current = tm_indices[chain_id][-1][1]
        length_current = tm_indices[chain_id][-1][2]

        # check before anything else to save computation time
        if desired_length <= length_current:

            # If there is no elongation to do, we skip to the next segment
            # and the coordinates of the ith segment are not modified
            return 0

        # (i-1)th TM Segment
        end_previous = tm_indices[chain_id][-2][1]

        # Compute the number of residues that are required to elongate the current segment
        elongation_left_to_do = desired_length - length_current



        # Randomly choose the number of residues to elongate downstream ( toward the C-terminal )
        downstream = random.randint(0, elongation_left_to_do)

        lefts = None

        # The new end of this final tm should not exceed the protein length
        if downstream + end_current > protein_length:

            # Hence take everyting that is between the end of the current tm segment and the start of the next tm segment
            new_end_coordinates = protein_length

            # What is " left " from downstream that could not be taken because the protein is too short after the last tm ? 
            lefts = downstream - (protein_length - end_current)


        else:

            new_end_coordinates = end_current + downstream


        upstream = elongation_left_to_do - downstream
        if lefts:

            upstream += lefts


        # The new start of this tm segment should not be lower than the end of the previous tm segment
        if start_current - upstream < end_previous:

            new_start_coordinates = end_previous + 1 

        else:

            new_start_coordinates = start_current - upstream     


        tm_indices[chain_id][-1] =(new_start_coordinates, new_end_coordinates, new_end_coordinates - new_start_coordinates + 1)

    return 0