function filtered_value = firstOrderFilter(alpha, new_value, previous_filtered_value)
    filtered_value = alpha * new_value + (1 - alpha) * previous_filtered_value;
end