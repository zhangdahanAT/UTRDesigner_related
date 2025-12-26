# -*- coding: utf-8 -*-
"""
Created on Fri Dec 26 08:15:59 2025

@author: Administrator
"""


import os

def split_file(file_path, part_size=100 * 1024 * 1024):
    """
    Splits a large file into smaller parts.

    Args:
    - file_path (str): The path to the file to be split.
    - part_size (int): The maximum size of each part in bytes. Default is 100MB.
    """
    # Get the base filename and extension
    base_name = os.path.basename(file_path)
    file_name, ext = os.path.splitext(base_name)
    
    # Open the file to read
    with open(file_path, 'rb') as f:
        part_number = 1
        while chunk := f.read(part_size):
            # Create a new part file with the format filename_partX.ext
            part_file_path = f"{file_name}_part{part_number}{ext}"
            with open(part_file_path, 'wb') as part_file:
                part_file.write(chunk)
            print(f"Created {part_file_path}")
            part_number += 1

def main():
    folder_path = r"D:\Backup\Documents\UTRdesigner\Arabidopsis DMS analysis"  # Update folder path if needed
    for file_name in os.listdir(folder_path):
        file_path = os.path.join(folder_path, file_name)
        if os.path.isfile(file_path):
            # Only process files larger than 100MB
            if os.path.getsize(file_path) > 100 * 1024 * 900:
                print(f"Splitting file: {file_name}")
                split_file(file_path)
            else:
                print(f"Skipping small file: {file_name}")
        else:
            print(f"Skipping directory: {file_name}")

if __name__ == "__main__":
    main()
