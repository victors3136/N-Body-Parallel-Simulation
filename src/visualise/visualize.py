import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import glob
import pandas as pd
import os


class NBodyVisualizer:
    def __init__(self, data_dir='csv'):
        if not os.path.exists(data_dir):
            os.makedirs(data_dir)

        self.files = sorted(glob.glob(f'{data_dir}/csv/positions_iter_*.csv'),
                            key=lambda x: int(x.split('_')[-1].split('.')[0]))
        if not self.files:
            raise FileNotFoundError(f"No position files found in {data_dir}! Run the simulation first.")

        print(f"Found {len(self.files)} position files")

        self.initial_data = pd.read_csv(self.files[0])
        self.n_particles = len(self.initial_data)
        print(f"Number of particles: {self.n_particles}")

        plt.style.use('dark_background')
        self.fig, self.ax = plt.subplots(figsize=(12, 12))

        self.scatter = self.ax.scatter(
            self.initial_data['x'],
            self.initial_data['y'],
            c='white',
            alpha=0.6,
            s=20  # Size of points
        )

        padding = 0.1e6  # 10% padding
        self.ax.set_xlim(-padding, 1e6 + padding)
        self.ax.set_ylim(-padding, 1e6 + padding)

        self.ax.set_title('N-Body Simulation', fontsize=14, pad=10)
        self.ax.set_xlabel('X Position')
        self.ax.set_ylabel('Y Position')

        self.ax.grid(True, alpha=0.2)

        self.iteration_text = self.ax.text(
            0.02, 0.98, '',
            transform=self.ax.transAxes,
            color='white',
            fontsize=12
        )

        print("\nInitial positions (first 5 particles):")
        print(self.initial_data.head())

    def update(self, frame):
        try:
            data = pd.read_csv(self.files[frame])

            if frame % 50 == 0:
                print(f"\nFrame {frame}:")
                print("Position ranges:")
                print(f"X: {data['x'].min():.2f} to {data['x'].max():.2f}")
                print(f"Y: {data['y'].min():.2f} to {data['y'].max():.2f}")

            self.scatter.set_offsets(data[['x', 'y']].values)

            self.iteration_text.set_text(f'Iteration: {frame}')

            return self.scatter, self.iteration_text

        except Exception as e:
            print(f"\nError in frame {frame}: {e}")
            print(f"File: {self.files[frame]}")
            return self.scatter, self.iteration_text

    def animate(self, interval=50, save_animation=False):
        print("\nStarting animation...")
        anim = animation.FuncAnimation(
            self.fig,
            self.update,
            frames=len(self.files),
            interval=interval,
            blit=True
        )

        if save_animation:
            print("\nSaving animation...")
            anim.save('nbody_simulation.gif',
                      writer='pillow',
                      fps=30)
            print("Animation saved as 'nbody_simulation.gif'")

        print("\nDisplaying animation...")
        plt.show()


if __name__ == "__main__":
    try:
        visualizer = NBodyVisualizer('../../outputs/threaded/')
        visualizer.animate(interval=20, save_animation=True)
    except Exception as e:
        print(f"\nCritical Error: {e}")

        import traceback

        traceback.print_exc()
