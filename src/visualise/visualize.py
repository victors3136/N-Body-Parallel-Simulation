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

        self.mass_scale = 500.0 / self.initial_data['mass'].max()

        plt.style.use('dark_background')
        self.fig, self.ax = plt.subplots(figsize=(12, 12))

        self.scatter = self.ax.scatter(
            self.initial_data['x'],
            self.initial_data['y'],
            s=self.initial_data['mass'] * self.mass_scale,
            c=self.initial_data['mass'],
            cmap='plasma',
            alpha=0.6,
            vmin=self.initial_data['mass'].min(),
            vmax=self.initial_data['mass'].max()
        )

        padding = 0.1e6
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

        self.mass_text = self.ax.text(
            0.02, 0.94, f'Mass Range: {self.initial_data["mass"].min():.0f} - {self.initial_data["mass"].max():.0f}',
            transform=self.ax.transAxes,
            color='white',
            fontsize=10
        )

        self.colorbar = self.fig.colorbar(self.scatter)
        self.colorbar.set_label('Mass', color='white')
        self.colorbar.ax.yaxis.set_tick_params(color='white')
        plt.setp(plt.getp(self.colorbar.ax.axes, 'yticklabels'), color='white')

        print("\nInitial positions and masses (first 5 particles):")
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
            self.scatter.set_sizes(data['mass'] * self.mass_scale)
            self.scatter.set_array(data['mass'])

            self.iteration_text.set_text(f'Iteration: {frame}')
            self.mass_text.set_text(f'Mass Range: {data["mass"].min():.0f} - {data["mass"].max():.0f}')

            return self.scatter, self.iteration_text, self.mass_text

        except Exception as e:
            print(f"\nError in frame {frame}: {e}")
            print(f"File: {self.files[frame]}")
            return self.scatter, self.iteration_text, self.mass_text

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
        visualizer = NBodyVisualizer('../../outputs/distributed/')
        visualizer.animate(interval=20, save_animation=True)
    except Exception as e:
        print(f"\nCritical Error: {e}")
        import traceback
        traceback.print_exc()
