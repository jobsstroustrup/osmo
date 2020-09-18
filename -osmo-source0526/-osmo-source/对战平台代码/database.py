#!/usr/bin/env python3

#####################################################
#                                                   #
#     ______        _______..___  ___.   ______     #
#    /  __  \      /       ||   \/   |  /  __  \    #
#   |  |  |  |    |   (----`|  \  /  | |  |  |  |   #
#   |  |  |  |     \   \    |  |\/|  | |  |  |  |   #
#   |  `--'  | .----)   |   |  |  |  | |  `--'  |   #
#    \______/  |_______/    |__|  |__|  \______/    #
#                                                   #
#                                                   #
#####################################################

# This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

import sqlite3
import time

class Database():
    def __init__(self, arg = None):
        self.connect = sqlite3.connect("data/" + str(round(time.time() * 1000)) + ".db")
        self.cursor = self.connect.cursor()

    def save_frame(self, frame_count, allcells):
        """Save the state of the cells in a frame.

        Args:
            frame_count: the index of current frame.
            allcells: list of all cells.
        Returns:
            

        """
        self.cursor.execute("""CREATE TABLE FRAME_{}
        (
            ID INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
            X  REAL    NOT NULL,
            Y  REAL    NOT NULL,
            VX REAL    NOT NULL,
            VY REAL    NOT NULL,
            R  REAL    NOT NULL
        );""".format(frame_count))
        for cell in allcells:
            if not cell.dead:
                self.cursor.execute(
                    "INSERT INTO FRAME_{} (X, Y, VX, VY, R) VALUES (?, ?, ?, ?, ?);".format(frame_count),
                    [cell.pos[0], cell.pos[1], cell.veloc[0], cell.veloc[1], cell.radius]
                )
        self.connect.commit()

    def save_game(self, data):
        """Save the state of the cells in the game.

        Args:
            data: data.
        Returns:
            

        """
        for i in range(len(data)):
            self.save_frame(i, data[i])
        self.connect.close()
